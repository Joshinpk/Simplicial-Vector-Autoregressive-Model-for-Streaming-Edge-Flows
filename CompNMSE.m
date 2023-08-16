
function nmse=CompNMSE(f_time_series,f_pred_time_series)
T=size(f_time_series,2);
E=size(f_time_series,1);
dif_sq=(f_time_series-f_pred_time_series).^2;
orig_sq=f_time_series.^2;

for t=1:T
    for i=1:E
        nmse(i,t)=sum(dif_sq(i,1:t))/sum(orig_sq(i,1:t));
    end
end
end