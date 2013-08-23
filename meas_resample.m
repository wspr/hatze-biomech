function calc = meas_resample(orig_meas,l,Norig,Nmeas,Ncalc)

orig_pos = linspace(0,l,Norig);
meas_pos = linspace(0,l,Nmeas);
calc_pos = linspace(0,l,Ncalc);

meas = interp1(orig_pos,orig_meas,meas_pos,'cubic');
calc = interp1(meas_pos,meas,calc_pos,'cubic');

end