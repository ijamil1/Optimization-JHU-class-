function [ofv,optsoln] = simplex1(A,b,c,basis)

szA = size(A);
rows = szA(1);
cols = szA(2);


pre_tab = rand(rows+1, cols+2);
for i = 1:(rows+1)
    for j = 1:(cols+2)
        if (i==1)
            if (j==1)
                preTab(i,j)=1;
            elseif (j < cols+2)
                preTab(i,j) = -1 * c(j-1);
            else
                preTab(i,j) = 0;   
            end
        else
            if (j == 1)
                preTab(i,j) = 0;
            elseif (j < cols+2)
                preTab(i,j) = A(i-1,j-1);
            else
                preTab(i,j) = b(i-1);
            end
        end
    end
end
B = A(:,basis); %assuming basisCols is horizontal vec/array
cb_t = [];
for i = 1:length(basis)
    cb_t = [cb_t basis(i)+1];
end
B_Inv = inv(B);
bottomPortion = B_Inv*preTab(2:rows+1,:);
preTab = [preTab(1,:);bottomPortion];
topRow = preTab(1,:) + ((-1*preTab(1,cb_t))*preTab(2:rows+1,:));
preTab = [topRow;preTab(2:rows+1,:)];

rn_t = preTab(1, 2:cols+1) > 0;
iteration = 1;
while (sum(rn_t)>0)
    rats(preTab)
    in_b=-1;
    max = 0;
    for i = 2:cols+1
        if preTab(1,i)>max
            max = preTab(1,i);
            in_b = i;
        end
    end
   
  %find which basis gets kicked out 
  minRatio = 1000000000;
  out_bIdx=0;
  for i = 2:rows+1
      if (preTab(i, in_b)<=0)
      else
          tempRatio = preTab(i,cols+2)/preTab(i,in_b);
          
          if (tempRatio < minRatio)
              minRatio = tempRatio;
              out_bIdx = i;
          end  
      end
  end
  
  preTab(out_bIdx,:) = preTab(out_bIdx,:)/preTab(out_bIdx,in_b);
  
  
  for i = 1:rows+1
      if (i~=out_bIdx)
          preTab(i, :) = preTab(i,:) + (-preTab(i,in_b) * preTab(out_bIdx,:));
      end
  end
  
  
  
  
  
%outBasisCol = basis(out_bIdx-1);
 % for j = 1:length(basis)
%     if (basis(j) == outBasisCol)
%     else
%         newbasis = [newbasis basis(j)+1];
%     end 
% end
% B = preTab(2:rows+1,newbasis)
% preTab = [preTab(1,:);(B\preTab(2:rows+1,:))];
% cb_t = newbasis;
% topRow = preTab(1,:) + ((-1*preTab(1,cb_t))*preTab(2:rows+1,:));
% newMat = [preTab(2:rows+1,:)];
% preTab = [topRow;newMat]

iteration = iteration + 1;
rn_t = preTab(1, 2:cols+1) > 0;
end
rats(preTab)
ofv = preTab(1,cols+2);
optsoln = preTab(2:rows+1,cols+2);

end
