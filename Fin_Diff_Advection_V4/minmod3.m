function m = minmod3(a,b,c)
T = a.*b > 0 & b.*c > 0;
m = zeros(size(a));
m(T) = sign(a(T)) .* min(min(abs(a(T)), abs(b(T))),abs(c(T)));
end