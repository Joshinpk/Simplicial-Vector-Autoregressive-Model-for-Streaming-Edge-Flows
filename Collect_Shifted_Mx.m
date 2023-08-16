function [Fs]=Collect_Shifted_Mx(f_val,L,MaxHops,tag)

scale=1;
fs=f_val;
if tag=='lwr'
    Fs=1*f_val;
    for k=1:MaxHops
        fs=scale*L*fs;
        Fs=[Fs,fs];
    end
elseif tag=='upr'
        Fs=[];
    for k=1:MaxHops
        fs=scale*L*fs;
        Fs=[Fs,fs];
    end
else
    dips('Error: invalid tag for Laplacian shift')
end
end