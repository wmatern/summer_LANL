function m = minmod(a,b)
T = a.*b > 0;
m = zeros(size(a));
m(T) = sign(a(T)) .* min( abs(a(T)), abs(b(T)) );
end