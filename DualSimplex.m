function [ofv,optvec] = DualSimplex(A,b,c,basis)

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
basic_vars = preTab(2:rows+1,cols+2); % actual basic vars of basic solution not nec feasible
bv = basic_vars < 0;
iteration = 1;
while (sum(bv)> 0)
    rats(preTab)
    in_b=-1;
    min = 1;
    row = -1;
    
    for i = 1:length(basic_vars)
        if basic_vars(i) < min
            min = basic_vars(i);
            row = i+1;
        end
    end
    
   % pivot 
   
  minRatio = 1000000000;
  pivCol=0;
  for i = 2:cols+1
      if (preTab(row, i)>=0)
      else
          tempRatio = -preTab(1,i)/-preTab(row,i);
          
          if (tempRatio < minRatio)
              minRatio = tempRatio;
              pivCol = i;
          end  
      end
  end
  
  preTab(row,:) = preTab(row,:)/preTab(row,pivCol);
  
  
  for i = 1:rows+1
      if (i~=row)
          preTab(i, :) = preTab(i,:) + (-preTab(i,pivCol) * preTab(row,:));
      end
  end
  
  
  
  
  

iteration = iteration + 1;
basic_vars = preTab(2:rows+1,cols+2); % actual basic vars of basic solution not nec feasible
bv = basic_vars < 0;
end
optvec = [];
for j = 1:cols
    if (preTab(1,j+1) ~= 0)
        optvec = [optvec 0];
    else
        for i = 1:rows
            if (preTab(i+1, j+1) == 1)
                optvec=[optvec preTab(i+1, cols+2)];
            end
        end
    
    end
end
rats(preTab)
optvec=transpose(optvec);
ofv = preTab(1,cols+2);


end