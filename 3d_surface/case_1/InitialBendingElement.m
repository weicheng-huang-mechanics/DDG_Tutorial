function bElement = InitialBendingElement(rodParams, triangle)
% This function initializes the bending elements for a discrete plate
%
%   Input:
%       rodParams - the defined plate struct contains the physical and
%                   numerical parameters of the simulated system
%       triangle  - the triangle mesh (nt * 3)
%
%   Output:
%       bElement  - Struct array of bending elements with 
%                   geometric and physical properties

EI = rodParams.EI;
x = rodParams.x;
nt = rodParams.nt;

bElement = struct('triIndex', {}, 'nodeIndex', {}, 'globalIndex', {}, 'nBar', {}, 'EI_local', {});

temp = 1;
for i = 1:nt
    for j = i+1:nt
        localT_1 = triangle(i,:);
        localT_2 = triangle(j,:);
        
        % check if connect
        connectNode = 0;
        connectPair = zeros(2,1);
        for ii = 1:3
            for jj = 1:3
                if ( localT_1(ii) == localT_2(jj)  )
                    connectPair(connectNode+1) = localT_1(ii);
                    connectNode = connectNode + 1;
                end
            end
        end
        
        if (connectNode == 2)
            bElement(temp).triIndex(1) = i;
            bElement(temp).triIndex(2) = j;
            
            bElement(temp).nodeIndex(2) = connectPair(1);
            bElement(temp).nodeIndex(3) = connectPair(2);
            
            for ii = 1:3
                if ( localT_1(ii) ~= connectPair(1) )
                    if ( localT_1(ii) ~= connectPair(2) )
                        bElement(temp).nodeIndex(1) = localT_1(ii);
                    end
                end
            end
            
            for ii = 1:3
                if ( localT_2(ii) ~= connectPair(1) )
                    if ( localT_2(ii) ~= connectPair(2) )
                        bElement(temp).nodeIndex(4) = localT_2(ii);
                    end
                end
            end
            
            % Compute Kappa Bar
            x_1 = getVertex(x, bElement(temp).nodeIndex(1));
            x_2 = getVertex(x, bElement(temp).nodeIndex(2));
            x_3 = getVertex(x, bElement(temp).nodeIndex(3));
            x_4 = getVertex(x, bElement(temp).nodeIndex(4));
            
            e_1 = x_3 - x_1;
            e_2 = x_2 - x_1;
            e_3 = x_2 - x_4;
            e_4 = x_3 - x_4;

            n_1 = cross(e_1, e_2);
            n_2 = cross(e_3, e_4);

            norm_1 = n_1 / norm(n_1);
            norm_2 = n_2 / norm(n_2);

            bElement(temp).nBar = norm(norm_1 - norm_2);
                        
            bElement(temp).globalIndex = zeros(11,1);
    
            bElement(temp).globalIndex(1) = 3 * (bElement(temp).nodeIndex(1)-1) + 1;
            bElement(temp).globalIndex(2) = 3 * (bElement(temp).nodeIndex(1)-1) + 2;
            bElement(temp).globalIndex(3) = 3 * (bElement(temp).nodeIndex(1)-1) + 3;
            
            bElement(temp).globalIndex(4) = 3 * (bElement(temp).nodeIndex(2)-1) + 1;
            bElement(temp).globalIndex(5) = 3 * (bElement(temp).nodeIndex(2)-1) + 2;
            bElement(temp).globalIndex(6) = 3 * (bElement(temp).nodeIndex(2)-1) + 3;
    
            bElement(temp).globalIndex(7) = 3 * (bElement(temp).nodeIndex(3)-1) + 1;
            bElement(temp).globalIndex(8) = 3 * (bElement(temp).nodeIndex(3)-1) + 2;
            bElement(temp).globalIndex(9) = 3 * (bElement(temp).nodeIndex(3)-1) + 3;
            
            bElement(temp).globalIndex(10) = 3 * (bElement(temp).nodeIndex(4)-1) + 1;
            bElement(temp).globalIndex(11) = 3 * (bElement(temp).nodeIndex(4)-1) + 2;
            bElement(temp).globalIndex(12) = 3 * (bElement(temp).nodeIndex(4)-1) + 3;
            
            bElement(temp).EI_local = EI;
            
            temp = temp + 1;
        end
        
    end
    
end

end
