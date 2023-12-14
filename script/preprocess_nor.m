function new_pupil = preprocess_nor(baseline,pupil,samplerate)
new_pupil = pupil - mean(pupil(:,(baseline(1)*samplerate/1000):(baseline(2)*samplerate/1000)),2);


end