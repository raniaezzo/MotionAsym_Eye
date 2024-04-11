# MotionAsym_Eye

INSTRUCTIONS:

- Change the setup.json file to reflect the correct paths of the data and the edf2asc files.
- Edf2asc conversion is proprietary and can be downloaded through the SR Research Portal.


~~~~~~~~~~~~~~~~~~~~~~~~~

NOTES/CHECK:

- add_blink.m and add_outside.m also have 500 and 1800 (stimulus period hard coded). Make sure that these values are not in terms of sampling rate and are instead milliseconds.
- Why are screen parameters loaded in for each file? Are they not the same across the board? (e.g., scr.disp_sizeX)

RECENT FIXES:
- convert.m had 1300 and 500 as hard coded values which go into the Dat_stim and Dat_all. I incorporated sampling rate.
- check_outside.m was hard coded. Fixed to make it dynamic (based on screen center) and DVA threshold of 1.5 for "outside" labelling.
- add_outside.m fixed the stimStart / stimEnd to be relative to sampling rate.
- add_blink.m fixed the stimStart / stipend in the same way.
- In convert, I am making Dat_all and Dat stim now add each trial to the end rather than beginning of the matrix. This might have been messing up some indexing of the Dat_all and Dat_stim matrices?
- I added edits to trial to downsample/upsample if needed.


~~~~~~~~~~~~~~~~~~~~~~~~~

preprocess_MS_01.m

For all subjects, for each motion direction/location, this script will call the  following scripts, then filters the timeseries (MATLAB's filtfilt.m) identifies Microsaccades (scripts below) and create a summary of Microsaccade properties. Inside the loop, multiple MS segments are joined into one matrix (within a trial) and saccades (defined by velocity and amplitude are turned to NaNs). It outputs (_events.mat) and (_rate.mat) files and corresponding .pngs. In the main script 8 files are saved: 

	convert.m: will convert EDF file to separate MAT files with timeseries (_DAT_all, _DAT_stim files), and timepoint labels (_tab). The timepoints are defined based on the printed MESSAGES inside the experimental code (e.g, TRIAL_START, TRIAL_END). _DAT_all contains all trial data (TRIAL_START to TRIAL_END per row) and leaves out data in-between. _DAT_stim contains all trial data from (STIM_ON to STIM OFF). It appears like this:
MSG	2333363 EVENT_ClearScreen
MSG	2333824 TRIAL_START 573
MSG	2333484 TRIAL_END

	Note about Eyelink bug: EDF2ASC conversion occasionally prints the TRIAL_END after the next TRIAL_START. When this happens, it occurs immediately after. If you try to re-arrange the MSG or ASC file based on timestamp order this does not fix the issue b/c sometimes it is received after. So a fix, messes up other trials. Best fix: Keep as is, and for the "broken trials" use the next TRIAL_START. I verified this is not due to fixation breaks etc.

	converge.m: will take the (_tab) file from above and add columns for trial behavior.

	check_outside.m: will take the Data_Stim period and check for eye position deviation from fixation

	add_outside.m: creates new (_tab) file which includes the boolean column above if any values exceeded the threshold from screen center during stimulus presentation.

	add_blink.m: updates file (_tab) to include a column with 0 or 1 indicating whether a blink occurred during the stimulus period based on Eyelink detection.

	findSamplingRate.m: returns the sampling rate according to that EDF/MSG file.
		NOTE: check that the outside and blink parameters account for any differences in sampling rate.

	omitblinks.m: turns a column with replicated timestamps and turns any blink period to nans. Note that the data is only during trial, so there are more blinks that reflected in the Dat_all_blink.mat files.

	segmentnonBlinks2.m: creates separate segments of the eyetrace between blinks. This will make a column in the trial matrix with 1s, 2s, 3s, etc corresponding to each eye trace between blinks.

	computevelocity.m: velocity computation is done separately for x,y. This is done AFTER the following steps for X and Y separately: (1) X and Y are converted from pixels to deg of visual angle. (2) fir1(35, 0.05) is used to create a low pass filter: This has a filter order of 35 (moderate # of filter coefficients), and 0.05 as the normalized cutoff frequency (5% of the Nyquist frequency; i.e is sampling rate = 1000 that would be 25 Hz. So frequencies below 25 Hz would be passed through. This is a type 1 linear-phase finite impulse response filter. (3) After this filtfilt(filter,1,x) performs a zero-phase digital filtering on X and Y (meaning the filter is applied forward and reverse directions to eliminate phase distortion). computevelocity.m effectively then computes the mean velocity for each point using a window of size 3 (considering neighboring points).

	microsaccMerge.m: detects monocular candidates of microsaccades based on relative velocity factor of (6x greater than median?), and minimum Microsaccade duration of 6ms, and merges micro saccades that are within 15ms of one another.  

	flterAM.m: filters the amplitude to be within 0.05 and 1 deg of visual angle.

	getms.m: filter MS matrix to leave out microsaccades that overlap with the end of the nSampleCutOff of 2500 ms.

	calulcate_direction.m: calculates direction of each microsaccades within stimulus period using 2 methods (components, and amplitude).
		calulcate_direction_RE.m (re-written for cleaner code and non-negative values. Inconsequential).

	count_ms.m: produces a summary table with all microsaccade characteristics for the session


~~~~~~~~~~~~~~~~~~~~~~~~~

preprocess_MS_01.m FILES

Tab files (_tab.mat contains all columns):
Col 1: trial attempt # (this is 1:N, including trials that are aborted). These numbers are not trialIDs, which would filter out mistrials.
Col 2: trial start based on print TRIAL_START (in eyelink sampling units)
Col 3: Unused
Col 4: trial sync time, usually the same as trial start? (Unused)
Col 5: stimulus start based on STIMULUS_ON (missing for some early subjects - but since stimulus was always 1300ms after TRIAL_START, can be retrieved)
Col 6: when EVENT_CLEAR_SCREEN prints
Col 7: Broke fixation based on 1.5 degree criteria in experiment
Col 8: TRIAL_END
Col 9: polar angle location
Col 10: motion direction
Col 11: tilt magnitude
Col 12: direction with tilt offset
Col 13: RT - 500
Col 14: Correct/Incorrect
Col 15: # of samples outside of the 1.5 distance from fixation
Col 16: # of samples during blink

Blink files (_blink.mat: this info is included in the tab file above)
Col 1: Start of the blink (in eyelink sampling units)
Col 2: End of the blink
Col 3: # of samples during blink (col 2 - col 1)

Outside files (_outside.mat: this info includes a boolean for each timestamp when eye position is 1.5 deg away from center)

Dat files (dat_stim for stimulus period only, dat_all for trial period only):
Col 1: Time sample
Col 2: X
Col 3: Y
Col 4: Pupil (arbitrary units)

Events File (_events.mat)
Row = Trial
Column = Binary value indicating whether Microsaccade or Not

Rate File (_rate.mat)
1 x 2500 file with rate

DIR_microsaccadeMatrix.mat (MS_TEMP) includes the Microsaccade data filtered to be within the correct amplitude range
*Also excludes microsaccades that intersect with the end of the nSampleCutOff which is set to 2500 ms after trial onset
Col 1: onset of saccade
Col 2: end of saccade
Col 3: peak velocity of saccade (vpeak)
Col 4: horizontal component     (dx)  differences in horizontal and vertical positions between the onset and offset points of each saccade
Col 5: vertical component       (dy)
Col 6: horizontal amplitude     (dX) differences in extreme horizontal and vertical positions between the onset and offset points of each saccade. Insights into the extent or range of movement
Col 7: vertical amplitude       (dY)
Col 8: microsaccade start (temporal onset) - redundant with Col 1
Col 9: microsaccade end (temporal onset) - redundant with Col 2
Col 10: trial ID number
Col 11: Amplitude (pythagorean distance based on Col 6 and 7)
Col 12: Direction from amplitude (col 6-7)
Col 13: Direction from components (col 4-5)

DIR_summary.mat (contains summary of microsaccade characteristics for a given session)
Rows = [# of MS in condition, mean AMP, max AMP, min AMP, mean duration, max duration, min duration, mean velocity peak, max velocity peak, min velocity peak]
Col 1: loc 315 (obl loc)
Col 2: loc 135 (obl loc)
Col 3: loc 225 (obl loc)
Col 4: loc 45 (obl loc)
Col 5: loc 270 (card loc)
Col 6: loc 90 (card loc)
Col 7: loc 180 (card loc)
Col 8: loc 0 (card loc)
Col 9: all microsaccades in session (sum)
*4 columns within a session should be 0s, since only 4 locations are tested per session.


~~~~~~~~~~~~~~~~~~~~~~~~~

plot_timeseriesPerCondition_02.m: plot the MS rate per condition per subject as specified in the script (either by tilt, direction, location, etc.). Also plots the median timeseries. Shaded error is error of the across session mean.
	*Used to be finalprocess16_RE.m, finalprocess17_RE.m, finalprocess17.m (replaced all)

	computeMSRateRT: creates a rate timeseries from the tab file and the MS file. Does so based on array of direction, location, and tilt values. Also computes overall RT (filtered by 2 s, and unfiltered RT. 
		* Used to be countigure_new.m

	shadedErrorBar.m: plots the standard error of the mean.


~~~~~~~~~~~~~~~~~~~~~~~~~~

compare_timeseries_03.m: plots the latency (defined as the time between last MS and the next relative to stimulus period), and plot them per condition.

	compute_latency() and fit_Line() functions defined within the script.

~~~~~~~~~~~~~~~~~~~~~~~~~~

new_script_for_Dat_all.m: 














~~~~~~~~~~~~~~~~~~~~~~~~~~

~ DO NEXT:

- Note: S03 VU for non_radial - there is block 1 and block 2 eyelink files. (TO DO)
		Block 1: expRes: 339 rows, but not every number is available. Goes from 20-22 (with 0 in the middle)Removed_Duplicate_
	S04 for radial - there are 2 blocks for every condition (TO DO)

finalprocess14.m 

new_script_for_Dat_all.m: makes Dat_all for pupil (trial start --> 2200)


~~~~~~~~~~~~~~~~~~~~~~~~~

sw_session.m

	Calls sw_withinsession or sw_withinsession2 or withsession: plot accuracy and MS rate

slope_figure.m

	Rate over time?

script.m : computes accuracy and reaction time with and without micro saccades.

plotcorr.m : computes correlation between MS rate and accuracy. 

Plot1123.m : plots polar histograms for each direction at each location.

plot_figure3.m : plots bars for # of MS per tilt angle, accuracy
	
	Calls acc_rec_count.m : computes tab for w and w/o MS for clockwise and counterclockwise

plot_figure2 : plots # of MS per tilt (deg_stacked)

plot_figure : plots # of MS per tilt (deg_stacked) - seems similar to above





Other sub functions: split_condition, 
	shadedErrorBar.m : probably for permutation?


~~~~~~~~~~~~~~~~~~~~~~~~~

REQUIRED SOFTWARE:

Matlab (Checked with R2021b)

SR Research Eyelink Developers Kit (must be installed on local computer)
	
