function [new_d, new_y] = dfunc(F, d, Y, h)

new_d = 0.99*d;

if rand < 0.5
    new_y = Y + 1;
else
    new_y = Y + 2;
end

if F > 50
    new_d = 1;
end

end