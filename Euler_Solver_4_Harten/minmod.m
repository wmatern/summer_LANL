function m = minmod(a,b)

m = 0.5 * (sign(a) + sign(b)).*min(abs(a),abs(b));

end