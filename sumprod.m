function [Z] = sumprod( A, B, w, its )

[hNodesnum, gNodesnum] = size(w);
msgMatrix = ones(hNodesnum, hNodesnum, gNodesnum);
oldMsgMatrix = msgMatrix;
beliefXi = ones(hNodesnum, gNodesnum);
beliefXiXj = zeros(hNodesnum, hNodesnum, gNodesnum, gNodesnum);


for it = 1: its
 for i= 1:hNodesnum %%for each node in H
    for j = 1:hNodesnum %% calc msg passed to all ne
        if (A(i,j) == 1)
            for xj= 1:gNodesnum %% calc msgs to each node in G
                   
                msgMatrix(i,j,xj) = getpx(i, j, xj, A, B, w, oldMsgMatrix);
             
            end
        end
    end
 end
 disp(msgMatrix)
 oldMsgMatrix = msgMatrix;  %% update old mas matrix
 %%calc beliefs xi
 for i = 1: hNodesnum
     for xi = 1:gNodesnum
         beliefXi(i, xi)  = exp(w(i,xi)) * getMsgPassBelf(i, xi, msgMatrix);
     end
 end

gamma = getNormGamma(beliefXi)
for i = 1: hNodesnum

        beliefXi(i,:) = beliefXi(i,:) / gamma(i);

end
% disp(beliefXi)
 %%calc beliefs for xi xj
 for i = 1: hNodesnum
     for j = 1: hNodesnum 
         if (A(i,j) == 1)
             for xi = 1: gNodesnum
                 for xj = 1: gNodesnum
                     beliefXiXj(i, j, xi, xj) = (1*B(xi,xj)) * exp(w(i, xi)) * exp(w(j, xj)) * getMsgPassedByNe(i, j, xi, msgMatrix, A) * getMsgPassedByNe(j, i, xj, msgMatrix, A); 
                 end
             end
         end
     end 
 end
end
disp(beliefXiXj);
for i = 1: hNodesnum

        beliefXiXj(i,:) = beliefXiXj(i,:) / gamma(i);

end
% disp(beliefXiXj);
%% calc Bethe energy
h1=0;
for i= 1: hNodesnum
    for xi = 1: gNodesnum
        h1 = h1 + beliefXi(i, xi)* log(beliefXi(i, xi));
    end
end
h2 =0;
b3 =0;
count =0;
for i = 1: hNodesnum
    for j = 1: i
        if (A(i,j) == 1)
            for xi = 1: gNodesnum
                for xj = 1 : gNodesnum
                    b1 = beliefXi(i, xi)*beliefXi(j, xj);
                    b2 = beliefXiXj(i, j, xi, xj)/b1;
                    if (b2 ~= 0)
                        b3 = log(b2);
                    end
                    h2 = h2 + beliefXiXj(i, j, xi, xj) * b3;
                    count = count +1;
                end
            end
        end
    end
end
%% calc H
disp (h1);
disp(h2);
%H = -h1-h2;
%disp(H);
u1 =0;
for i= 1: hNodesnum
    for xi = 1: gNodesnum
        if (exp(w(i, xi)) == 0)
            c1 = 1;
        else
            c1 = log(exp(w(i, xi)));
        end
        u1 = u1 + beliefXi(i, xi)* c1;
    end
end
u2=0;
for i = 1: hNodesnum
    for j = 1: hNodesnum
        for xi = 1: gNodesnum
            for xj = 1 : gNodesnum
                if (B(xi,xj) == 0)
                    b2 = 1;
                else
                    b2 = log(1*B(xi,xj));   %%check %%
                end
                u2 = u2 + beliefXiXj(i, j, xi, xj) * b2;
            end
        end
    end
end
%%calc U
disp(u1)
disp(u2)
%U = -u1-u2;
%disp(U);
%disp(U)
%logZ = U -  H
logZ = u1 + u2  - h1 - h2;
Z = exp(logZ)
end

 function prod = getpx(i, j, xj, A, B, w, oldMsgMatrix)
 [rowsA, colsA] = size(B);
 temp =0;
  for xi = 1:colsA 
    temp = temp + exp(w(i, xi)) * 1*(B(xi, xj)) * getMsgPassedByNe(i, j, xi, oldMsgMatrix, A);
  end
  prod = temp;
 end
 
 function multiMsg = getMsgPassedByNe(i, j, xi, oldMsgMatrix, A)
    m =1;
    [rowsA, colsA] = size(A);
        for col = 1: colsA  %%find edges connected to i except j
            if (col ~= j) && (A(i, col) ==1)
            %%if (A(i, col) ==1)
                m = m * oldMsgMatrix(i, col, xi);
            end
        end
        multiMsg = m;
 end
 
 function m = getMsgPassBelf(i, xi, msgMatrix)
 [d1, d2, d3] = size(msgMatrix);
 m =1;
    for x = 1: d1
            m = m * msgMatrix(x, i, xi);    
    end
 end
 
 function p = getpxNew(i, j, xi, xj, A, B, w, msgMatrix)
    p = exp(w(i, xi)) * 1*(B(xi, xj)) * getMsgPassedByNe(i, j, xi, oldMsgMatrix, A);
 end
 
 function n = getNormFact(i, hNodesnum, msgMatrix, A)
 sumMatN = sum (msgMatrix, 3);
    for j = 1: hNodesnum
        if(A(i,j) == 1)
            n(i) = sumMatN(i, j);
        end
    end 
 end
 function gamma = getNormGamma(beliefXi)
    gamma = sum(beliefXi, 2)
 end
 
