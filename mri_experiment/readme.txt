ISSS_multiband pilot study
% % % % %
Changelog:
V1: Stimuli presentation, scanner synchrony
  6/21 Created directories. Started work on v1. Ported from ISSS_test. 
  # Finished adapting LoadStimuli.m, WaitForScannerTrigger.m, instructions.txt, DisplayInstructions_bkfw.m
  # Decided to encode subject data as follows:
    001_MH_isss_ ... 
  # File outputs now have unique names based on subj #, initials, and scan protocol
  # Subjects can proceed forwards and backwards through instructions without an immediate error on incorrect
    button press during first instruction. 
  6/22 Experiment can now run!
  # Fixed error in LoadStimuli that would not output sound correctly (stimuli are mono, output is stereo. 
    Renamed LoadStimuli to LoadStimuli_mono to show the difference. 
  # Checked that output and timing are still working (they are). 
  # Overhauled LoadStimuli_mono so it generates the answer, event, and jitter keys (LoadStimuli_mono_anskey)
  # Created answer key (1 = female, 0 = male)
  # Added fixation cross
  6/26 Finished V1
  # Changed pathing and order of WaitForScannerTrigger function, hoping to increase timing accuracy of stimuli. 
  6/29 Tested script in scanner, made minor changes
V2: Collecting subject response
  6/26 Started V2
  # Added results folder to collect all outputs, organized by subject number
  6/28 Copied over to my flash drive, pushed to GitHub, and copied direc to Box to ensure preparedness for Thursday.
  6/29 Copied small changes from test, and added a few features:
    Fixed screen and headphone ID problems
    Cursor now hides after screen opens
    Cross is smaller
    Made sure Skip Instructions was working properly
    Added "waiting for experimenter" screen when waiting for first pulse
    Made fixation cross to appear after first pulse
  # Added try statement to catch errors easier
  # Experimenters can now press escape to end code!
V3: Careful constraints on stimuli now in place:
  - 16 specific sentences from a list of 32 [now 64 - 7/3] (counterbalanced across runs?)
  - 8 noise
  - 8 stimuli
  6/29 Put silent and 1ch stimuli into stimuli directory
  # Rewrote Loadstimuli_keys so now it loads ALL stimuli and creates a key that plays the above stimuli
  # Removed ALL references to mulitple runs for this experiment, most variables are no longer cells, just matrices. 
  6/30 Added RTBox scripts to isss_multiband study folder
  # Added a shortcut to open RTBox directory
  # Added option to simulate keyboard as RTBox if RTBox isn't connected
  # Rewrote DisplayInstructions to run from RTBox commands
  7/3 Changed some pathing of old scripts, unused stimuli, to decrease file size
  # Adjusted LoadStimuli to load stimuli correctly. 
  # Renamed LoadStimuli to LoadStimuliAndKeys. Now counterbalances stimuli and outputs stimuli duration. 
  7/5 Double-checked stimuli, integrated stimulicheck to experiment script
  # Locked timing to the end of stimuli presentation. 
  # Changed timing variable names, fixed issues with timing. Finished V3?
  # Timing was broke. Will fix in v4? 
V4: Completed 7/6 testing, many bugs. Will fix in V4:
  # Fix timing errors (DONE 7/7)
  # Check hide cursor
  # Bigger fixation cross
V5: New, shorter blocks
  # Blocks will now contain only 16 total stimuli:
    8 sentences (chosen from 8 sentence structures)
    4 silent
    4 noise
  # Twice as many blocks, each of half length
  # Found loud stimuli to test headphone volume
  # Data now outputs in separate function, as txt and xls
  # Added training block (~2 minutes, ten unique stimuli) 
  # Make instructions that make sense (go over with Lee?)
  # Now using equalized stimuli
  # Added resting state/DTI blocks (black screen, no cross, ppt slide)
  # Introduced changes from 7/11 test, updated scan notes
  # Practice block runs on laptop 
RDT V1: Started coding RDT
  # Same 4.000 second silent window, stimuli are 1.5 seconds long with 0.2 second gap between
  # Jitter is now 0.800 seconds maximum. 
  # Stimuli need to be counterbalanced at 50% same, 50% different. 
  # Practice will be 5 trials long. 
  # 16 trials consisting of pairs of stimuli, two runs, means 64 stimuli?
  # LoadStimuliAndKeys_rdt should be done.
  # New stimuli have been made (RMS complete) and subsequently discarded, still need stimuli
  # Jitter is back, complete with variable length stimuli. 
  # Results now output correctly. 
First subject: Errors in Rhythm code associated with timing (same as we had earlier with isss_multiband)
  $ Stop scripts from quitting between runs.
  $ Update instructions to include random press on noise trials. 
  $ Add training component (lang task) to use within the scan to help assess volume of stimuli
      ~6 trials, with multiband imaging
      Have feedback
Subject 3: Scripts are completely ready to go. 
  
  
Future versions: 
  $ Preprocessing results
    $ How to do STC for multiband? 
  $ Data analysis! 