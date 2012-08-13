function m = minmod3(a,b,c)

for i=1:length(a(:,1))
    for j=1:length(a(1,:))
        
        if a(i,j)*b(i,j) > 0 && b(i,j)*c(i,j) > 0
            m(i,j) = sign(a(i,j)) * min( abs([a(i,j) b(i,j) c(i,j)]) );
        else
            m(i,j) = 0.0;
        end
    end
end

end