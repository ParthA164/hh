cost=[2 3 11 7; 1 7 6 1; 5 8 15 9];
A=[14 16 5];
B=[6 10 15 4];

if sum(A)==sum(B)
    fprintf("Balanced");
else
    fprintf("unbalanced");
    if sum(A)<sum(B)
        cost(end+1,:)=zeros(1,size(B,2));
        A(end+1)=sum(B)-sum(A);
    elseif sum(B)<sum(A)
        cost(:,end+1)=zeros(1,size(A,2));
        B(end+1)=sum(A)-sum(B);
    end
end

%%
icost=cost;
X=zeros(size(cost));
[m,n]=size(cost);
for i=1:m
    for j=1:n
        hh=min(cost(:))
        [rowind,colind]=find(hh==cost)
        x11=min(A(rowind),B(colind))
        [val,ind]=max(x11)
        ii=rowind(ind)
        jj=colind(ind)

        y11=min(A(ii),B(jj))
        X(ii,jj)=y11
        A(ii)=A(ii)-y11
        B(jj)=B(jj)-y11

        cost(ii,jj)=Inf
    end
end

%%printing transportation table
transportation=array2table(X)
disp(transportation)


%%printing least cost
value=sum(sum(icost.*X))
fprintf("Least Cost is %d",value)




clc
clear all
a=[2 1 1 0 0;2 5 0 1 0; 2 3 0 0 1];
b=[50; 100; 90];
cost=[4 10 0 0 0 0];
A=[a b];
Var={'x1', 'x2', 's1', 's2', 's3','sol'};
%%defining basi`c variables and calculating zjcj
bv=[3 4 5];
zjcj=cost(bv)*A-cost
%%displaying initial simplex table
simplex_table=[zjcj;A]
array2table(simplex_table,'VariableNames',Var)

RUN=true;
while RUN
if any(zjcj(1:end-1)<0) %%check for negative value
    fprintf('The current BFS is not optimal \n')
    zc=zjcj(1:end-1)
    [Enter_val, pvt_col]=min(zc)


if all(A(:,pvt_col)<=0)
    error('LPP is unbounded'); %%for ratio denominator can't be 0;
else
    sol=A(:,end)
    column=A(:,pvt_col)
    for i=1:size(A,1)
        if column(i)>0
            ratio(i)=sol(i)./column(i)
        else 
            ratio(i)=inf
        end
    end
    [leaving_val, pvt_row]=min(ratio)
end

%%entering and exiting variable
bv(pvt_row)=pvt_col   %%replaced leaving variable with entering variable
pvt_key=A(pvt_row,pvt_col)
A(pvt_row,:)=A(pvt_row,:)./pvt_key
%A(pvt_row,:)=A(pvt_row,:)./pvt_key;

%row operation
for i=1:size(A,1)
    if i~=pvt_row
        A(i,:)=A(i,:)-A(i,pvt_col).*A(pvt_row,:);
    end
end
zjcj=cost(bv)*A-cost;
next_table=[zjcj;A];
array2table(next_table,'VariableNames',Var)
else
    RUN=false;
    fprintf('The table is optimal\n');
    z=input('Enter 0 for minimization and 1 for max \n'); 
        if z==0
          Obj_value=-zjcj(end)
        else
            Obj_value=zjcj(end)
        end
        fprintf('The final optimal value is %f\n', Obj_value);
end
end



clc
clear all
a=[2 1 1 0 0;2 5 0 1 0; 2 3 0 0 1];
b=[50; 100; 90];
cost=[4 10 0 0 0 0];
A=[a b];
Var={'x1', 'x2', 's1', 's2', 's3','sol'};
%%defining basi`c variables and calculating zjcj
bv=[3 4 5];
zjcj=cost(bv)*A-cost
%%displaying initial simplex table
simplex_table=[zjcj;A]
array2table(simplex_table,'VariableNames',Var)

RUN=true;
while RUN
if any(zjcj(1:end-1)<0) %%check for negative value
    fprintf('The current BFS is not optimal \n')
    zc=zjcj(1:end-1)
    [Enter_val, pvt_col]=min(zc)


if all(A(:,pvt_col)<=0)
    error('LPP is unbounded'); %%for ratio denominator can't be 0;
else
    sol=A(:,end)
    column=A(:,pvt_col)
    for i=1:size(A,1)
        if column(i)>0
            ratio(i)=sol(i)./column(i)
        else 
            ratio(i)=inf
        end
    end
    [leaving_val, pvt_row]=min(ratio)
end

%%entering and exiting variable
bv(pvt_row)=pvt_col   %%replaced leaving variable with entering variable
pvt_key=A(pvt_row,pvt_col)
A(pvt_row,:)=A(pvt_row,:)./pvt_key
%A(pvt_row,:)=A(pvt_row,:)./pvt_key;

%row operation
for i=1:size(A,1)
    if i~=pvt_row
        A(i,:)=A(i,:)-A(i,pvt_col).*A(pvt_row,:);
    end
end
zjcj=cost(bv)*A-cost;
next_table=[zjcj;A];
array2table(next_table,'VariableNames',Var)
else
    RUN=false;
    fprintf('The table is optimal\n');
    z=input('Enter 0 for minimization and 1 for max \n'); 
        if z==0
          Obj_value=-zjcj(end)
        else
            Obj_value=zjcj(end)
        end
        fprintf('The final optimal value is %f\n', Obj_value);
end
end



clc
clear all
a=[-1 -1 1 1 0;-1 2 -4 0 1];
b=[-5; -8];
cost=[-2 0 -1 0 0 0];
A=[a b];
Var={'x1', 'x2', 'x3', 's1', 's2', 'sol'};
%%defining basi`c variables and calculating zjcj
bv=[4 5];
zjcj=cost(bv)*A-cost
%%displaying initial simplex table
simplex_table=[zjcj;A]
array2table(simplex_table,'VariableNames',Var)

Run=true
while Run
    sol=A(:,end);
    ratio=[];
    if any(sol<0)
        fprintf('The current BFS is not feasable\n');
        [leaving_variable,pvt_row]=min(sol);
        fprintf('Leaving row %d \n', pvt_row);
        for i=1:size(A,2)-1
            if A(pvt_row,i)<0
                ratio(i)=abs(zjcj(i)/A(pvt_row,i));
            else
                ratio(i)=inf;
            end
        end
        [entering_variable, pvt_col]=max(-ratio);
        fprintf('Entering column %d \n', pvt_col)
        bv(pvt_row)=pvt_col   %%replaced leaving variable with entering variable
        pvt_key=A(pvt_row,pvt_col)
        A(pvt_row,:)=A(pvt_row,:)./pvt_key
        %row operation
        for i=1:size(A,1)
            if i~=pvt_row
                A(i,:)=A(i,:)-A(i,pvt_col).*A(pvt_row,:);
            end
        end
        zjcj=cost(bv)*A-cost;
        next_table=[zjcj;A];
        array2table(next_table,'VariableNames',Var)
        else
            RUN=false;
            fprintf('The table is optimal\n');
            z=input('Enter 0 for minimization and 1 for max \n'); 
            if z==0
                Obj_value=-zjcj(end)
            else
                Obj_value=zjcj(end)
            end
            fprintf('The final optimal value is %f\n', Obj_value);
    end
end

   


syms x1 x2 
%to declare x1 and x2 as symbols , not numbers or arrays

% define objective function
f1=x1-x2 + 2*x1^2 + 2*x1*x2 + x2^2;
fx=inline(f1);   %to convert to function check fx(1,1)
fobj=@(x) fx(x(:,1),x(:,2));

% Gradient
grad=gradient(f1);
G=inline(grad);
gradx=@(x) G(x(:,1),x(:,2));

% Hessian
h=hessian(f1);
hx=inline(h);

% initialization
x0=[1 1];
maxiter=4; % maximum iterations
tot=1e-3; % set the tolerance to 0.001
iter=0;   %initial iteration
X=[];

% Iterative steepest descent

while norm(double(gradx(x0))) > tot && iter < maxiter
    X=[X,x0];
    S=-double(gradx(x0));    % make sure gradient is numeric
    hval=double(hx(x0));     %make sure hessian is numeric
    lam=(S' * S) / (S' * hval * S)  % step size
    x0= x0 + lam * S'   %update step
    iter= iter + 1;
end

fprintf("Optimal solution x=[%f %f] \n", x0(1), x0(2));
fprintf("Optimal Value f(x)=%f \n",fobj(x0));
