function cutoff = CHopping(coeffs, tol) 
% cutoff has two number
% for coeffs = (c_0 , c_2 , ... ,c_j , ... , c_j2 , ... , c_n)
% (c_j , ... , c_j2) is the plateau , then cutoff = [j ; j2])


if ( tol >= 1 ) % input not enough
    cutoff = [1 ; 1];
    return
end

n = length(coeffs);
cutoff = [2 * n ; 2 * n]; %assume there is no plateau
if ( n < 17 )  %length of coeffs is too short to build a plateau
    return
end

b = abs(coeffs);
m = b(end)*ones(n, 1);
for j = n-1:-1:1
    m(j) = max(b(j), m(j+1));
end   
if ( m(1) == 0 )
    cutoff = [1 ; 1];
    return
end
envelope = m/m(1);

for j = 2:n
    j2 = round(1.25*j + 5); 
    if ( j2 > n )
        % there is no plateau: exit
        return
    end      
    e1 = envelope(j);
    e2 = envelope(j2);
    r = 3*(1 - log(e1)/log(tol));
    plateau = (e1 == 0) | (e2/e1 > r);
    if ( plateau )
        cutoff = [j ; j2];
        break
    end
end
end