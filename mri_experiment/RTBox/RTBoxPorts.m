function [out, vers] = RTBoxPorts(allPorts)
% [availPorts, vers]= RTBoxPorts(); % return available RTBox ports and vers
% [allPorts, vers]= RTBoxPorts(1); % all RTBox ports and vers (close if needed)
% 
% [port, st]= RTBoxPorts(busyPorts); 
% cell input If input is a cellstr, it will be treated as in-use RTBox ports.
% Then this function will open 1st available RTBox port, and return the port
% name as 1st output, and the 2nd output will be a struct containing serial
% handle, RTBox version, clock unit, latency timer and host MAC address. In case
% no available RTBox is found, the port will be empty, and the 2nd struct output
% will show available and busy ports.

% 090901 Wrote it (XL)
% 100102 use non-block write to avoid problem of some ports
% 120101 add second output
% 120701 take care of return from bootloader or RTBboxADC
% 170426 RTBox.m, RTBoxClass.m etc start to call this to open port
% 170502 Include MACAddress() and LatencyTimer() here for convenience
% 170508 filters FTDI ports for all OS.

if nargin<1, allPorts = false; end
toOpen = ischar(allPorts) || iscellstr(allPorts); % to open first avail port 

verbo = IOPort('Verbosity', 0); % shut up screen output and error
cln = onCleanup(@() IOPort('Verbosity', verbo));

if ~toOpen && logical(allPorts) % ask all ports, so close
    try RTBox('CloseAll'); end
    try RTBoxClass.instances('closeAll'); end
    try RTBoxSimple('close'); end
end

if     ispc,   ports = ftdi_vcp_win();
elseif ismac,  ports = ftdi_vcp_mac();
elseif isunix, ports = ftdi_vcp_lnx();
else, error('Unsupported system: %s.', computer);
end

out = {}; vers = [];
cfgStr = 'BaudRate=115200 ReceiveTimeout=0.2 PollLatency=0';
rec = struct('avail', '', 'busy', ''); % for error message only

for i = 1:numel(ports)
    port = ports{i};
    if any(strcmp(port, allPorts)), continue; end % avoid multi-open in unix
    s = IOPort('OpenSerialPort', port, cfgStr);
    if s<0, rec.busy{end+1} = port; continue; end
    
    idn = RTBox_idn(s);
    if strncmp(idn, '?', 1) % maybe in boot/RTBoxADC
        IOPort('Write', s, 'R', 0); % return to application
        WaitSecs('YieldSecs', 0.1); drawnow;
        idn = RTBox_idn(s);
    end
    if numel(idn) < 21 % re-open to fix rare ID failure
        IOPort('Close', s); WaitSecs('YieldSecs', 0.1); drawnow;
        s = IOPort('OpenSerialPort', port, cfgStr);
        idn = RTBox_idn(s);
        if numel(idn)<21, idn = RTBox_idn(s); end % try one more time
    end
    if numel(idn)==21 && strncmp(idn, 'USTCRTBOX', 9)
        v = str2double(idn(19:21));
        if v>100, v = v/100; end % v510, rather than v5.1
        if toOpen
            out = port; break;            
        else
            out{end+1} = port; vers(end+1) = v; %#ok<*AGROW>
        end
    else
        rec.avail{end+1} = port; % avail but not RTBox
    end
    IOPort('Close', s); % close it
end

if isempty(out), vers = rec; return; end % no RTBox found
if ~toOpen, return; end

% The rest is for RTBox.m, RTBoxClass.m etc to set up serial port
try [oldVal, err] = LatencyTimer(port, 2); lat = min(2, oldVal);
catch me, oldVal = 16; err = me.message; lat = 16; % in case of error
end
if ~isempty(err) % error, failed to change, or change not effective
    if oldVal>2
        warning('LatencyTimer:Fail', ['%s\nThis simply means failure to speed ' ...
            'up USB-serial port reading. It won''t affect RTBox function.'], err);
    end
    lat = oldVal;
end
if oldVal>lat % close/re-open to make change effect
    IOPort('Close', s);
    s = IOPort('OpenSerialPort', port, cfgStr);
end

vers = struct('ser', s, 'version', v, 'clockUnit', 1/str2double(idn(11:16)));
vers.latencyTimer = lat/1000;
vers.MAC = [0 MACAddress];

function idn = RTBox_idn(s) % return RTBox idn str
IOPort('Purge', s); % clear buffer
IOPort('Write', s, 'X', 0); % non-blocking to avoid problem for some ports
idn = char(IOPort('Read', s, 1, 21)); % USTCRTBOX,921600,v6.1

function ports = ftdi_vcp_lnx()
% Return FTDI serial ports under Linux.
ports = dir('/dev/ttyUSB*');
ports = {ports.name};
for i = numel(ports):-1:1
    ports{i} = ['/dev/' ports{i}];
    [err, str] = system(['udevadm info -a -n ' ports{i}]);
    if err, continue; end % give up for now
    ind = regexp(str, '{idProduct}=+"6001"', 'once');
    if isempty(ind), ports(i) = []; continue; end
    i1 = regexp(str, char([10 10]));
    i2 = i1(find(i1>ind, 1)); i1 = i1(find(i1<ind, 1, 'last'));
    ind = regexp(str(i1:i2), '{idVendor}=+"0403"', 'once');
    if isempty(ind), ports(i) = []; end
end

function ports = ftdi_vcp_mac()
% Return FTDI serial ports under OSX.
ports = dir('/dev/cu.usbserial*');
ports = {ports.name};
[~, str] = system('system_profiler SPUSBDataType');
hex4 = '[0-9a-fA-F]{4}';
for i = numel(ports):-1:1
    ports{i} = ['/dev/' ports{i}];
    p = ports{1}(end-7:end);
    p = p([7:8 5:6 3:4 1:2]); % Location ID
    ind = regexp(str, ['Location ID:\s*0x' p], 'ignorecase');
    if numel(ind) ~= 1, continue; end % give up for now
    i1 = regexp(str, char([10 10]));
    i2 = i1(find(i1>ind, 1)); i1 = i1(find(i1<ind, 1, 'last'));
    ip = regexp(str(i1:i2), ['(?<=Product ID:\s*0x)' hex4], 'match', 'once');
    iv = regexp(str(i1:i2), ['(?<=Vendor ID:\s*0x)' hex4], 'match', 'once');
    if ~strcmp(ip, '6001') || ~strcmp(iv, '0403'), ports(i) = []; end
end

function ports = ftdi_vcp_win()
% Return FTDI serial ports under Windows.
HLM = 'HKEY_LOCAL_MACHINE';
sub = 'HARDWARE\DEVICEMAP\SERIALCOMM';
try % winqueryreg is fast: ~1 ms
    ports = {};
    nams = winqueryreg('name', HLM, sub); % list \Device\serial&VCP
    for i = 1:numel(nams)
        if isempty(strfind(nams{i}, 'VCP')), continue; end
        ports{end+1} = winqueryreg(HLM, sub, nams{i});
    end
catch % if no winqueryreg, fallback to reg query for Octave: ~60 ms
    [~, str] = system(['reg.exe query ' HLM '\' sub]);
    ports = regexp(str, 'VCP\d*\s+REG_SZ\s+(COM\d{1,3})', 'tokens');
    for i = 1:numel(ports), ports{i} = ports{i}{1}; end
end

ftdi = 'HKEY_LOCAL_MACHINE\SYSTEM\CurrentControlSet\Enum\FTDIBUS';
[err, str] = system(['reg.exe query ' ftdi ' /s /f PortName']);
if err % early reg.exe
    [~, str] = system(['reg.exe query ' ftdi ' /s | findstr PortName']);
end
ftdiPorts = regexp(str, 'COM\d{1,3}(?=\s)', 'match');
for i = numel(ports):-1:1
    if ~any(strcmp(ports{i}, ftdiPorts)), ports(i) = []; continue; end
    ports{i} = ['\\.\' ports{i}];
end

function val = regquery(key, name)
% Query the registry value of 'name' in 'key'.
% Call winqueryreg if available, and use reg query otherwise.
try
    i = strfind(key, '\'); i = i(1);
    val = winqueryreg(key(1:i-1), key(i+1:end), name);
    if isnumeric(val), val = double(val); end
catch
    [~, str] = system(['reg.exe query "' key '" /v ' name]);
    tok = regexp(str, [name '\s+(REG_\w+)\s+(\w+)'], 'tokens', 'once');
    if numel(tok)<2, val = []; return; end
    val = tok{2};
    if ~isempty(strfind(tok{1}, '_SZ')), return; end % char type
    if strncmp(val, '0x', 2) % maybe always hex type?
        val = sscanf(val, '0x%x');
    else
        val = str2double(val);
    end
end

function mac = MACAddress()
% Return computer MAC address in uint8 of length 6. If all attemps fail, there
% will be a warning, and last 6 char of host name will be returned for RTBox.
% sprintf('%02X-%02X-%02X-%02X-%02X-%02X', MACAddress) % '-' separated hex
mac = zeros(1, 6, 'uint8');
try %#ok<*TRYNC> OSX and Linux will return from this block
    a = Screen('computer');
    hex = a.MACAddress;
    a = sscanf(hex, '%2x%*c', 6);
    if numel(a)==6, mac(:) = a; return; end
end

try % java approach faster than getmac, mainly for Windows
    ni = java.net.NetworkInterface.getNetworkInterfaces;
    while ni.hasMoreElements
        a = ni.nextElement.getHardwareAddress;
        if numel(a)==6 && ~all(a==0)
            mac(:) = typecast(a, 'uint8'); % from int8
            return; % 1st is likely ethernet adaptor
        end
    end
end

try % system command is slow
    if ispc, cmd = 'getmac'; else, cmd = 'ifconfig'; end
    [~, str] = system(cmd);
    a = '[0-9a-fA-F]{2}'; expr = ['(' a '[:-]{1}){5}' a];
    hex = regexp(str, expr, 'match', 'once');
    a = sscanf(hex, '%2x%*c', 6);
    if numel(a)==6, mac(:) = a; return; end
end

warning('MACAddress:Fail', 'Using last 6 char of hostname as MACaddresss');
[~, nam] = system('hostname'); nam = strtrim(nam);
if numel(nam)<6, nam = ['myhost' nam]; end
mac(:) = nam(end+(-5:0));

function [val, errmsg] = LatencyTimer(port, msecs)
% Query/change FTDI USB-serial port latency timer. 
%  lat = LatencyTimer(port); % query only 
%  val = LatencyTimer(port, msecs); % query and set to msecs if val>msecs
% 
% Administrator/sudo privilege is normally needed to change the latency timer.
errmsg = '';
warnID = 'LatencyTimer:RestrictedUser';
warnmsg = 'Failed to change latency timer due to insufficient privilege.';
if ispc
    port = strrep(port, '\\.\', '');
    ftdi = 'HKEY_LOCAL_MACHINE\SYSTEM\CurrentControlSet\Enum\FTDIBUS';
    [err, str] = system(['reg.exe query ' ftdi ' /s /v PortName']);
    if err, [~, str] = system_file(['reg.exe query ' ftdi ' /s']); end % win XP?
    expr = ['HKEY_.*(?=\n\s+PortName\s+REG_SZ\s+' port ')'];
    key = regexp(str, expr, 'match', 'dotexceptnewline', 'once');
    val = regquery(key, 'LatencyTimer');
    if nargin<2, return; end
    msecs = uint8(msecs); % round it and make it within 255
    if val <= msecs, return; end
    
    fid = fopen('temp.reg', 'w'); % create a reg file
    fprintf(fid, 'REGEDIT4\n[%s]\n"LatencyTimer"=dword:%08x\n', key, msecs);
    fclose(fid);
    % change registry, which will fail if not administrator
    [err, txt] = system('reg.exe import temp.reg 2>&1');
    delete('temp.reg');
    if err
        errmsg = [warnmsg ' ' txt 'You need to start Matlab/Octave by right-clicking'...
            ' its shortcut or executable, and Run as administrator.'];
        if nargout<2, warning(warnID, WrapString(errmsg)); end
    end
elseif ismac
    % port not needed. After the change, all FTDI serial ports may be affected.
    useFTDI = true; % use FTDI driver
    folder = '/Library/Extensions/FTDIUSBSerialDriver.kext'; % for later OS
    fname = fullfile(folder, '/Contents/Info.plist');
    if ~exist(fname, 'file')
        folder = '/System/Library/Extensions/FTDIUSBSerialDriver.kext';
        fname = fullfile(folder, '/Contents/Info.plist');
    end
    if ~exist(fname, 'file')
        useFTDI = false; % use driver from Apple: different keys
        folder = '/System/Library/Extensions/IOUSBFamily.kext/Contents/PlugIns/AppleUSBFTDI.kext';
        fname = fullfile(folder, '/Contents/Info.plist');
    end
    if ~exist(fname, 'file')
         error('LatencyTimer:plist', 'Info.plist not found.');
    end
    
    fid = fopen(fname);
    str = fread(fid, '*char')';
    fclose(fid);
    
    if useFTDI
        ind = regexp(str, '<key>FTDI2XXB', 'once');
        if isempty(ind), ind = regexp(str, '<key>FT2XXB', 'once'); end
    else
        ind = regexp(str, '<key>AppleUSBEFTDI-6001', 'once');
    end
    if isempty(ind)
        error('LatencyTimer:key', 'Failed to detect FTDI key.');
    end
    i2 = regexp(str(ind:end), '</dict>', 'once') + ind; % end of ConfigData
    expr = '<key>LatencyTimer</key>\s+<integer>\d{1,3}</integer>';
    [mat, i0, i1] = regexp(str(ind:i2), expr, 'match', 'start', 'end', 'once');
    if isempty(i0)
        % TODO:
        % if isempty(i0) && ~useFTDI, Insert LatencyTimer key
        error('LatencyTimer:key', 'Failed to detect LatencyTimer key.');
    end
    valStr = regexp(mat, '\d{1,3}(?=</integer>)', 'match', 'once');
    val = str2double(valStr);
    if nargin<2, return; end % query only
    msecs = uint8(msecs);
    if val <= msecs, return; end
    
    tmp = strrep(fname, '/Info.plist', '/tmpfoo');
    fid = fopen(tmp, 'w+'); % test privilege
    if fid<0
        fprintf(' You will be asked for sudo password to change the latency timer.\n');
        fprintf(' Enter to skip the change.\n');
        err = system('sudo -v');
        if err
            errmsg = warnmsg;
            if nargout<2, warning(warnID, WrapString(errmsg)); end
            return;
        end
    else
        fclose(fid);
        delete(tmp);
    end
    
    i0 = i0+ind-1; i1 = i1+ind-1; % index of mat in str, including
    mat = strrep(mat, [valStr '</integer>'], [num2str(msecs) '</integer>']);
   
    tmp = '/tmp/tmpfoo';
    fid = fopen(tmp, 'w+');
    fprintf(fid, '%s', str(1:i0-1)); % before mat
    fprintf(fid, '%s', mat); % modified mat
    fprintf(fid, '%s', str(i1+1:end)); % after mat
    fclose(fid);
    system(['sudo mv -f ' tmp ' ' fname]);
    system(['sudo touch ' folder]);
    system('sudo -k');
    errmsg = 'The change will take effect after you reboot the computer.';
    if nargout<2, warning([mfile ':rebootNeeded'], errmsg); end
else % for linux, no vendor related info needed, only the port name
    port = strrep(port, '/dev/', '');
    str = sprintf('cd /sys/bus/usb-serial/devices/%s;', port);
    [err, lat] = system([str 'cat latency_timer']); % query
    if err % unlikely happen
        error('LatencyTimer:readFail', ['Failed to read latency timer: ' lat]);
    end
    val = str2double(lat);
    if nargin<2, return; end % query only
    msecs = uint8(msecs);
    if val <= msecs, return; end
    
    system([str 'echo ' num2str(msecs) ' > latency_timer']);
    [~, lat] = system([str 'cat latency_timer']); % check for sure
    if msecs ~= str2double(lat)
        errmsg = [warnmsg ' You need to run Matlab as superuser.'];
        if nargout<2, warning(warnID, WrapString(errmsg)); end
    end
end

function [err, out] = system_file(cmd)
% The same as [err, out] = system(cmd), but better performance with large out.
persistent fname deleteFile
if isempty(deleteFile)
    fname = [tempdir 'reg_output_junk.txt'];
    deleteFile = onCleanup(@() delete(fname)); % delete only when matlab exit
end
err = system([cmd ' >' fname]);
fid = fopen(fname);
out = fread(fid, '*char')';
fclose(fid);
