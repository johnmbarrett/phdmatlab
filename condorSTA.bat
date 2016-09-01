if %1 lss 7 (
	set /a channel = %1 + 11
) else (
	if %1 gtr 54 (
		set /a channel = %1 + 27
	) else (
		set /a x = %1+1
		set /a y = x/8
		set /a z = y+1
		set /a w = x%%8
		set /a channel = 10*z + w + 1
	)
)

set "experiment_folder=V:\retina\John B\Electrophysiolgy\Wild-Type Stim\%2"
set "data_folder=%experiment_folder%\%3_%~4"
set "spike_file=times_%~4_channel_%channel%_MCD_trimmed_spikes.mat"
set "timings_file=%~4_photodiode_timings.mat"
set "stim_file=%~4.mj2"

copy "%data_folder%\%spike_file%" .
copy "%data_folder%\%timings_file%" .
copy "%experiment_folder%\%stim_file%" .

condorSTA.exe "%spike_file%" "%timings_file%" "%stim_file%" %5 %6