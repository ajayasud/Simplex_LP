% Code:
%% Stage 1: To convert the Problem into standard form and also to add artificial variables if necessary PLUS REDUNDANCY REMOVED

% coefficient matrix in its raw form is given as input:

A = [

         0    3.0000    2.0000    2.0000   -1.0000         0         0         0         0         0    1.0000         0         0         0         0         0
         0    1.0000    0.5000    2.0000    1.0350   -1.0000         0         0         0         0   -1.0300    1.0000         0         0         0         0
         0    1.8000   -1.5000    1.8000         0    1.0350   -1.0000         0         0         0         0   -1.0300    1.0000         0         0         0
         0   -0.4000   -1.5000   -1.0000         0         0    1.0350   -1.0000         0         0         0         0   -1.0300    1.0000         0         0
         0   -1.8000   -1.5000   -1.0000         0         0         0    1.0350   -1.0000         0         0         0         0   -1.0300    1.0000         0
         0   -1.8000   -0.2000   -1.0000         0         0         0         0    1.0350   -1.0000         0         0         0         0   -1.0300    1.0000
    1.0000   -5.5000    1.0000   -6.0000         0         0         0         0         0    1.0350         0         0         0         0         0   -1.0300
         0         0         0         0    1.0000         0         0         0         0         0         0         0         0         0         0         0
         0         0         0         0         0    1.0000         0         0         0         0         0         0         0         0         0         0
         0         0         0         0         0         0    1.0000         0         0         0         0         0         0         0         0         0
         0         0         0         0         0         0         0    1.0000         0         0         0         0         0         0         0         0
         0         0         0         0         0         0         0         0    1.0000         0         0         0         0         0         0         0
         0         0         0         0         0         0         0         0         0    1.0000         0         0         0         0         0         0
         0    1.0000         0         0         0         0         0         0         0         0         0         0         0         0         0         0
         0         0    1.0000         0         0         0         0         0         0         0         0         0         0         0         0         0
         0         0         0    1.0000         0         0         0         0         0         0         0         0         0         0         0         0
]

v = [0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1]                            % A vector 'v' holding depicting the signs of the constraints;
                                                                                        % (-1 for >=), (0 for =) and (1 for <=)
[m n] = size(A);            % this will return # of rows to 'm' and # of columns to 'n'

mat = zeros(m,1);           % this will create a zero vector of dim m*1


for k=1:m
if v(k)== 0
mat(k) = 1;                     % If the sign is '=', it adds an artificial variable; this is done by adding column 
                                %with a an entry = 1 corresponding to that constraint
A = [A mat];
mat = zeros(m,1);
end
end


for k=1:m
if v(k)== -1
mat(k) = 1;                     % If the sign is '>=' it addds an artificial variable
A = [A mat];                    %with a an entry = 1 corresponding to that constraint
mat = zeros(m,1);
end
end

[ba ab] = size(A)		% to define cost vector for the first phase, we will be using 
                        %this as all the artificial variables required have been added.

for k=1:m
if v(k)== 1
mat(k) = 1;
A = [A mat];            % if the sign is '<=' it adds a column with an entry = 1 corresponding to that constraint.
mat = zeros(m,1);       %Thus a slack variable is added.
end
end

[bbb ccc] = size(A) 		%to define indices of columns in initial basis we wil be using this;
                            % Now, we have created an initial Basis which
                            % is identity

for k=1:m
if v(k)== -1                
mat(k) = -1;            % if the sign is '<=' it adds a column with an entry = -1 corresponding to that constraint.
A = [A mat];            %Thus a slack variable is added.
mat = zeros(m,1);
end
end


A;                          % We now have a matrix A that is ready to be used in phase 1

% Redundancy Check : checking condition included at the last
[row col] = size(A)
B = A
for i = 1 : row
for j = 1 : col
B(i,j) = Inf;                           % We are definig a matrix B of the dimensions same as A with entries as INFINITY
end
end
count_red = 0                           % We define an arbitrary variable cout_red as ZERO
for t = 1 : row
	for r = t+1 : row
	u = zeros(col,1);                    % We are also defining a zero vector which will hold the RATIOS of values 
                                        % present in two rows that are being checked 
	for p = 1:col
	u(p,1) = (A(r,p)/A(t,p)) ;           
	end
	k = u(1)
		for v = 1 : col	
		if u(v) == k
		count_red = count_red + 1;       % We are checking if all the ratios are same
		end
		end

		if count_red == col
			if k >=1
				for w = 1 : col
				B(r,w) = 0  ;            % If the ratio > 1, it means the row with higher index is redundant
                                        % If the ratios are same corresponding row in B is assigned ZERO
				B;
				count_red=0;
				end
			else
				for w = 1 : col
				B(t,w) = 0 ;             % If the ratio > 1, it means the row with lower index is redundant
                                        % If the ratios are same corresponding row in B is assigned ZERO
				B;
				count_red=0;
				end
			end
		end
		count_red = 0;
	end
end


% Checkig condition : if any particular row is zero in B matrix, corresponding row in A is redundant
z = zeros(row, 1)
for i = 1 : row
if B(i,:)==0            % We check if any of the row in B has all entries as ZERO; If so corresponding row in A is redundant
z(i) = i
end
end

for i = 1 : row
if z(i) > 0
	A(z(i),:) = []      % The row in A that corresponds to the row in B that has all zeros, is finally removed here.
	z(i) = 0
	z = z - 1
end
end
A
% End of redundancy check

[mm nn] = size(A);


b = [2;.5;.4;.38;.36;.34;.3;2;2;2;2;2;2;1;1;1];         % b vector and cost vector, defined by 'cc' are given as inputs

cc = [-1;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];

c = zeros(nn,1);        
c(n+1 : ab) = 1;

iB = [n+1 : ccc]        % We are defining the the columns that should enter the initial Basis

tt = iB(1) - 1;         % This will be used later to remove constraints from A that are redundant at the end of phase 1.
mB = A(:, iB);          % Initial Basis is extracted from A
[row_count col_count] = size(A);  

%% Stage 2
enter = 5;              % We assign enter to be the maximum value of reduced cost coefficient(denoted by rcc) at a later stage
                        % If enter >0, it goes to the next iteration; else
                        % it returns the optimal solution. This is
                        % implementd here with while. For the first
                        % iteration to be started, we give an arbitrary
                        % value to 'enter'
while enter > 0.000000001

x = zeros(col_count,1);
x(iB) = inv(mB)*b;          % BFS is calculated





rcc = c(iB)'*inv(mB)*A-c' ;    % rcc is reduced cost coeffiecient and it is calculated in this step

enter = max(rcc);              % enter is given max value present in rcc

if enter < 0.000000001
fprintf('current basis is optimal');    % Optimality check: It enter < 0, it retrns the optimal solution 
                                        % and optimal objective value(Zopt); Else the simplex continues
x
Zopt = c(iB)'*inv(mB)*b
Zopt
break
end



%% Stage 3: Checking if theproblem is degenerate
deg_check = x(iB);
[~,basis_variable_count] = size(iB);
count = 1;                              % We check if any of Basis variables is zero
for s = 1 : basis_variable_count        % We define an arbritary variable (count) as 1 before starting the degenracy check     
	if(deg_check(s) == 0)               % if there is any Basis variables that is zero, count will be incremented by 1
                                        % If count > 1 we say it's
                                        % degenerate.
	count = count + 1;
	end
end
count


% If the problem is degenerate, the code uses Bland's rule to select a variable to enter the basis as per the following:
if count > 1
for u = 1 : col_count
	if rcc(u) > 0.00000001
	entering_index = u;     % We are finding the first index of RCC that is greater than zero.
	enter = rcc(u)
	break
	end
end
else
for t = 1 : col_count
	if rcc(t)== enter
	entering_index = t;        % If not degenerate, normal simplex is happening here
	break
	end
end
end

% End of this degeneracy code and simplex continues



%% Stage 4 : We check whether the LP is bounded or not.

ubc = zeros(col_count,1);                          % ubc is unbounded check

ubc(iB) = inv(mB)*A(:,entering_index);

var = max(ubc);                 % We assign max value of ubc to var

if var < 0.000000001            % if var < 0, it's unbounded;
fprintf('LP is unbounded');
break
end


for e= 1:col_count
		if ubc(e) <= 0.0000001                 % If the Lp is bounded we need to find the ratio for all POSITIVE values of ubc 
                                               % and to choose the minimum
                                               % of that.. So assign all
                                               % NON POSITIVE values to
                                               % INFINITY in ubc
		ubc(e) = Inf;
		end
end
                                                % here we the ratio and
                                                % also the minimum of that.
for d = 1:col_count
if ubc(d) > 0.000000001 && ubc(d) ~= inf
ubc(d) = x(d)/ubc(d);
end
end

leave= min(ubc);

%% Stage 5 If the problem is degenerate, the code uses Bland's rule to select a variable to leavethe basis as per the following:
if count > 1                            
for r = 1 : col_count
	if ubc(r)== leave
                                        % if there is degeneracy leave = 0
	exiting_index= r;                   % We are finding the first index of ubc whose rato is zero is greater than zero.
	break
	end
end
for t=1:row_count
	if iB(t) == exiting_index
	k = iB(t);
	iB(t) = entering_index;             % Here the index of entering column replaces the index of leaving column in Basis 
	break
	end
end
else
[leave exiting_index] = min(ubc);
for t=1:row_count
	if iB(t) == exiting_index
	k = iB(t);                          % Else normal simplex is happening here and index of entering column 
                                        % replaces the index of leaving column in Basis 
	iB(t) = entering_index;
	break
	end
end
end

% End of this degeneracy code and simplex continues


iB;
mB = A(:, iB);

end


%% Stage 6: PHASE 1 to PHASE 2

if Zopt > 0
fprintf('The given LP is infeasible\n')   % At the end of Phase 1, we are checking if LP s feasible

else
% if the LP is feasible, we start the building of start of Phase 2
fprintf('The given LP is feasible\n')
count_iB = 0
z = zeros(col_count,1);
for j = 1 : col_count
if c(j) == 1                            
z(j) = j;                           
end
end
                                    % These two for loops serve ONE purpose. That is to remove the columns corresponding to 
                                    % to artificial variables in 'A'
for i = 1 : col_count
if z(i) > 0
	A(:,z(i)) = [];
	z(i) = 0;
	z = z - 1;
	count_iB = count_iB + 1;
end
end



v = zeros(col_count,1);
for j = 1 : col_count
if c(j) == 1
v(j) = j;
end
end
                                    % These two for loops serve ONE purpose. That is to remove the COLUMNS corresponding to 
                                    % to artificial variables in 'C'
for i = 1 : col_count
if v(i) > 0
	c(v(i),:) = [];
	v(i) = 0;
	v= v - 1;
end
end

A
b
c = cc



for t = 1: row_count
if iB(t) > tt && iB(t) <= count_iB + tt		% We are checking if any of the artificial variable is present in Optimal solution of PHASE 1
                                            % If so, corresponding row in A
                                            % is redundant and hence we are
                                            % removing the constraint as as
                                            % whole.
A(t,:) = [];
iB(t) = [];
b(t) = [];
break
end
end

[row_count col_count] = size(A);
for t = 1: row_count
if iB(t) > tt                               % In the above discussed scenario we have to remove the column corresponding to the 
                                            % artificial variable that is
                                            % present in Basis.
iB(t) = iB(t) - count_iB;
end
end
iB;
A;
b;
c;



%% Stage 7: Phase 2
% We are finding if the solution at the end of phase 1 after removing all artificial variables is optimal or not. IF not optimal, we proceed with 2nd phase.

mB = A(:, iB);
[row_count col_count] = size(A);            

x = zeros(col_count,1);
x(iB) = inv(mB)*b
rcc = c(iB)'*inv(mB)*A-c' 
%rcc is reduced cost coeffiecient

enter = max(rcc)

if enter < 0.000000001
fprintf('current basis is optimal');
x
Zopt = c(iB)'*inv(mB)*b
Zopt

else
while enter > 0.000000001

x = zeros(col_count,1);
x(iB) = inv(mB)*b;                      % BFS is calculated





rcc = c(iB)'*inv(mB)*A-c' ;             % rcc is reduced cost coeffiecient and it is calculated in this step

enter = max(rcc);                       % enter is given max value present in rcc

if enter < 0.000000001
fprintf('current basis is optimal');    % Optimality check: It enter < 0, it retrns the optimal solution 
                                        % and optimal objective value(Zopt); Else the simplex continues
x
Zopt = c(iB)'*inv(mB)*b
Zopt
fprintf('Therefore, Maximised WSDMs net worth is %d \n', (-1 *Zopt))
break
end



%% Stage 8: Checking if theproblem is degenerate
deg_check = x(iB);
[~,basis_variable_count] = size(iB);
count = 1;                              % We check if any of Basis variables is zero
for s = 1 : basis_variable_count        % We define an arbritary variable (count) as 1 before starting the degenracy check     
	if(deg_check(s) == 0)               % if there is any Basis variables that is zero, count will be incremented by 1
                                        % If count > 1 we say it's
                                        % degenerate.
	count = count + 1;
	end
end
count


% If the problem is degenerate, the code uses Bland's rule to select a variable to enter the basis as per the following:
if count > 1
for u = 1 : col_count
	if rcc(u) > 0.00000001
	entering_index = u;                         % We are finding the first index of RCC that is greater than zero.
	enter = rcc(u)
	break
	end
end
else
for t = 1 : col_count
	if rcc(t)== enter
	entering_index = t;                         % If not degenerate, normal simplex is happening here
	break
	end
end
end

% End of this degeneracy code and simplex continues



%% Stage 9 : We check whether the LP is bounded or not.

ubc = zeros(col_count,1);                      % ubc is unbounded check

ubc(iB) = inv(mB)*A(:,entering_index);

var = max(ubc);                                % We assign max value of ubc to var

if var < 0.000000001                           % if var < 0, it's unbounded;
fprintf('LP is unbounded');
break
end


for e= 1:col_count
		if ubc(e) <= 0.0000001                 % If the Lp is bounded we need to find the ratio for all POSITIVE values of ubc 
                                               % and to choose the minimum
                                               % of that.. So assign all
                                               % NON POSITIVE values to
                                               % INFINITY in ubc
		ubc(e) = Inf;
		end
end
                                                % here we the ratio and
                                                % also the minimum of that.
for d = 1:col_count
if ubc(d) > 0.000000001 && ubc(d) ~= inf
ubc(d) = x(d)/ubc(d);
end
end

leave= min(ubc);

%% Stage 10: If the problem is degenerate, the code uses Bland's rule to select a variable to leavethe basis as per the following:
if count > 1                            
for r = 1 : col_count
	if ubc(r)== leave
                                        % if there is degeneracy leave = 0
	exiting_index= r;                   % We are finding the first index of ubc whose rato is zero is greater than zero.
	break
	end
end
for t=1:row_count
	if iB(t) == exiting_index
	k = iB(t);
	iB(t) = entering_index;             % Here the index of entering column replaces the index of leaving column in Basis 
	break
	end
end
else
[leave exiting_index] = min(ubc);
for t=1:row_count
	if iB(t) == exiting_index
	k = iB(t);                          % Else normal simplex is happening here and index of entering column 
                                        % replaces the index of leaving column in Basis 
	iB(t) = entering_index;
	break
	end
end
end

% End of this degeneracy code and simplex continues


iB;
mB = A(:, iB);

end


end                                     % for the if just above second while

end                                     % for the if to move from phase 1 to phase 2

