function m = minmod(a,b,c)

if a*b > 0 && b*c > 0
   m = sign(a) * min( abs([a b c]) );
elsme
   m = 0.0;
end

end