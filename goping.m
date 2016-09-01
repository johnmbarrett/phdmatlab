function goping
    sinwave = @(x,f) sin(2*pi*x*f/8192);
    sound(0.5*sinwave(1:8192/4,440))
end