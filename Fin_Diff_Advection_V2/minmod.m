function m = minmod(a,b)

for i=1:length(a(:,1))
    for j=1:length(a(1,:))
        
        if a(i,j)*b(i,j) > 0 > 0
            m(i,j) = sign(a(i,j)) .* min( abs([a(i,j) b(i,j)]) );
        else
            m(i,j) = 0.0;
        end
    end
end

end