function bElement = InitialBendingElement(rodParams, bend, sElement)

x = rodParams.x;

% build bending element
bElement = struct('edgeIndex', {}, 'nodeIndex', {}, 'globalIndex', {}, 'nodePos_1', {}, 'nodePos_2', {}, 'nodePos_3', {}, ...
    'voroLen', {}, 'theta_1', {}, 'theta_2', {}, 't_1', {}, 'd_11', {}, 'd_12', {}, 't_2', {}, 'd_21', {}, 'd_22', {}, ... 
    't_1_old', {}, 'd_11_old', {}, 'd_12_old', {}, 't_2_old', {}, 'd_21_old', {}, 'd_22_old', {}, ...
    'm_11', {}, 'm_12', {}, 'm_21', {}, 'm_22', {}, 'refTwist', {}, 'refTwist_old', {}, ... 
    'directSign_1', {}, 'directSign_2', {}, 'kappaBar', {}, 'EI_local', {}, 'GJ_local', {});

for i = 1:rodParams.nb
    bElement(i).edgeIndex = bend(i,:);
    
    % edge used in bending
    localEdge_1 = sElement(bElement(i).edgeIndex(1));
    localEdge_2 = sElement(bElement(i).edgeIndex(2));
    
    % local node index
    index1 = localEdge_1.nodeIndex(1);
    index2 = localEdge_1.nodeIndex(2);
    index3 = localEdge_2.nodeIndex(1);
    index4 = localEdge_2.nodeIndex(2);
    
    if (index1 == index3)
        bElement(i).nodeIndex(1) = index2;
        bElement(i).nodeIndex(2) = index1;
        bElement(i).nodeIndex(3) = index4;
    end
    
    if (index1 == index4)
        bElement(i).nodeIndex(1) = index2;
        bElement(i).nodeIndex(2) = index1;
        bElement(i).nodeIndex(3) = index3;
    end
    
    if (index2 == index3)
        bElement(i).nodeIndex(1) = index1;
        bElement(i).nodeIndex(2) = index2;
        bElement(i).nodeIndex(3) = index4;
    end
    
    if (index2 == index4)
        bElement(i).nodeIndex(1) = index1;
        bElement(i).nodeIndex(2) = index2;
        bElement(i).nodeIndex(3) = index3;
    end
    
    % node position and theta
    bElement(i).nodePos_1 = getVertex(x, bElement(i).nodeIndex(1));
    bElement(i).nodePos_2 = getVertex(x, bElement(i).nodeIndex(2));
    bElement(i).nodePos_3 = getVertex(x, bElement(i).nodeIndex(3));
    
    bElement(i).theta_1 = 0.0;
    bElement(i).theta_2 = 0.0;
    
    % build voroniLength
    bElement(i).voroLen = (norm(bElement(i).nodePos_2 - bElement(i).nodePos_1) + norm(bElement(i).nodePos_3 - bElement(i).nodePos_2) ) / 2;
    
    % build frist frame
    bElement(i).t_1 = (bElement(i).nodePos_2 - bElement(i).nodePos_1) / norm(bElement(i).nodePos_2 - bElement(i).nodePos_1);
    t_temp = [0;0;1];
    d1Tmp = cross(bElement(i).t_1, t_temp);
    if (abs(d1Tmp) < 1.0e-6)
        t_temp = [0;1;0];
        d1Tmp = cross(bElement(i).t_1, t_temp);
    end
    bElement(i).d_11 = d1Tmp;
    bElement(i).d_11 = bElement(i).d_11 / norm(bElement(i).d_11);
    bElement(i).d_12 = cross(bElement(i).t_1, bElement(i).d_11);
    bElement(i).d_12 = bElement(i).d_12 / norm(bElement(i).d_12);
    
    bElement(i).t_1_old  = bElement(i).t_1;
    bElement(i).d_11_old = bElement(i).d_11;
    bElement(i).d_12_old = bElement(i).d_12; 
    
    % build second frame
    bElement(i).t_2 = (bElement(i).nodePos_3 - bElement(i).nodePos_2) / norm(bElement(i).nodePos_3 - bElement(i).nodePos_2);
    bElement(i).d_21 = parallel_transport(bElement(i).d_11, bElement(i).t_1, bElement(i).t_2);
    bElement(i).d_21 = bElement(i).d_21 / norm(bElement(i).d_21);
    bElement(i).d_22 = cross(bElement(i).t_2, bElement(i).d_21);
    bElement(i).d_22 = bElement(i).d_22 / norm(bElement(i).d_22);
    
    bElement(i).t_2_old  = bElement(i).t_2;
    bElement(i).d_21_old = bElement(i).d_21;
    bElement(i).d_22_old = bElement(i).d_22;
    
    % build material frame
    cs = cos( bElement(i).theta_1 );
    ss = sin( bElement(i).theta_1 );
    bElement(i).m_11 =   cs * bElement(i).d_11 + ss * bElement(i).d_12;
    bElement(i).m_12 = - ss * bElement(i).d_11 + cs * bElement(i).d_12;
    
    cs = cos( bElement(i).theta_2 );
    ss = sin( bElement(i).theta_2 );
    bElement(i).m_21 =   cs * bElement(i).d_21 + ss * bElement(i).d_22;
    bElement(i).m_22 = - ss * bElement(i).d_21 + cs * bElement(i).d_22;
    
    % compute undeformed curvature
    kb = 2 * cross(bElement(i).t_1, bElement(i).t_2) / (1 + dot(bElement(i).t_1, bElement(i).t_2));
    bElement(i).kappaBar(1) =   0.5 * dot(kb, bElement(i).m_12 + bElement(i).m_22 );
    bElement(i).kappaBar(2) = - 0.5 * dot(kb, bElement(i).m_11 + bElement(i).m_21 );
    
    % Compare the sign with edge element
    if (dot(localEdge_1.t, bElement(i).t_1) >= 0)
        bElement(i).directSign_1 = 1.0;
    else
        bElement(i).directSign_1 = -1.0;
    end
    
    if (dot(localEdge_2.t, bElement(i).t_2) >= 0)
        bElement(i).directSign_2 = 1.0;
    else
        bElement(i).directSign_2 = -1.0;
    end
    
    % set reference twist
    bElement(i).refTwist = 0.0;
    bElement(i).refTwist_old = 0.0;
    
    % set global index
    bElement(i).globalIndex = zeros(11,1);
    
    bElement(i).globalIndex(1) = 3 * (bElement(i).nodeIndex(1)-1) + 1;
    bElement(i).globalIndex(2) = 3 * (bElement(i).nodeIndex(1)-1) + 2;
    bElement(i).globalIndex(3) = 3 * (bElement(i).nodeIndex(1)-1) + 3;
    
    bElement(i).globalIndex(4) = 3 * rodParams.nv + bElement(i).edgeIndex(1);
    
    bElement(i).globalIndex(5) = 3 * (bElement(i).nodeIndex(2)-1) + 1;
    bElement(i).globalIndex(6) = 3 * (bElement(i).nodeIndex(2)-1) + 2;
    bElement(i).globalIndex(7) = 3 * (bElement(i).nodeIndex(2)-1) + 3;
    
    bElement(i).globalIndex(8) = 3 * rodParams.nv + bElement(i).edgeIndex(2);
    
    bElement(i).globalIndex(9)  = 3 * (bElement(i).nodeIndex(3)-1) + 1;
    bElement(i).globalIndex(10) = 3 * (bElement(i).nodeIndex(3)-1) + 2;
    bElement(i).globalIndex(11) = 3 * (bElement(i).nodeIndex(3)-1) + 3;
    
    
    if ( norm(bElement(i).kappaBar) > 1e-2 )
        bElement(i).EI_local = 100 * rodParams.EI;
        bElement(i).GJ_local = 100 * rodParams.GJ;
    else
        bElement(i).EI_local = rodParams.EI;
        bElement(i).GJ_local = rodParams.GJ;
    end
end

end