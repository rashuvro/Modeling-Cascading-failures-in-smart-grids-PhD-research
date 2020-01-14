function new_h = hfunc(F, d, h)

new_h = 1.01 * h;

if F > 50
    new_h = 0;
end

end