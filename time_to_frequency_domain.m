function [X] = time_to_frequency_domain(x,dt,frequency_array,time_shift,complex_number)
number_of_time_steps = size(x,2);
number_of_frequencies = size(frequency_array,2);
X=zeros(1,number_of_frequencies);
w=2.*pi*frequency_array;
for n=1:number_of_time_steps
    t=n*dt+time_shift;
    X=X+x(n)*exp(-complex_number*w*t);
end
X=X*dt;
