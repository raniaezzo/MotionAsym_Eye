function MS_new=getms(MS,samplerate,nSampleCutOff)
    time = samplerate/1000;
    fil=MS(:,9)/time < nSampleCutOff;
    MS_new=MS(fil,:);
end