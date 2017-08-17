%===========================
% Delta by likelihood method
clear
r=0.05; % annual interest
sig=0.3; % volatility
K=35; % stike price
S_0=100; % initial stock price
T=4/12; %maturity in 4 months, which is 4/12 of the year
n= 88;  % there are approx. 88 trading days in 4 months
dt=T/n; % time step
% Simulate N sample paths of the stock process

K=20
result1=[];
for j=1:10

N=10^5; Delta=nan(N,1);
for i=1:N
    Z=randn(1,n);
    path=(r-sig^2/2)*dt+sig*sqrt(dt)*Z;
    path=cumprod([S_0,exp(path)]);
    A=mean(path); % average path/stock price
    jj=Z(1)/(S_0*sig*sqrt(dt));
    Delta(i)=exp(-r*T)*max(0,A-K)*jj;
end
D=mean(Delta);
stderr1=std(Delta)/sqrt(N);
result1=[result1; K D stderr1];
K=K+10
end
disp ('    K         Delta(LH)  stderr');
result1

% Delta by Pathwise Derivative method
clear
r=.05; % annual interest
sig=0.3; % volatility
K=35; % stike price
S_0=100; % initial stock price
T=4/12; %maturity in 4 months, which is 4/12 of the year
n= 88;  % there are approx. 88 trading days in 4 months
dt=T/n; % time step
% Simulate N sample paths of the stock process

result2=[];
K=20
for j=1:10

N=10^5; Delta=nan(N,1);Vega=Delta;
for i=1:N
    path=(r-sig^2/2)*dt+sig*sqrt(dt)*randn(1,n);
    path=cumprod([S_0,exp(path)]);
    A=mean(path); % average path/stock price
    Delta(i)=exp(-r*T)*(A>=K)*A/S_0;
end
D=mean(Delta);
stderr = std(Delta)/sqrt(N);
result2=[result2; K D stderr];
K=K+10
end
disp ('     K        Delta(LH)  stderr');
result2

% Delta by Finite Differences Method
clear
r=.05; % annual interest
sig=0.3; % volatility
K=35; % stike price
S_0=100; % initial stock price
T=4/12; %maturity in 4 months, which is 4/12 of the year
n= 88;  % there are approx. 88 trading days in 4 months
dt=T/n; % time step
% Simulate N sample paths of the stock process

h=1;  %to be used in formula (f(x+h) - f(x))/h
N=10000;
K=20
result3=[];
for K=20:10:110
    %find f(x)  in the terms (f(x+h) - f(x))/h
    ranvec= randn(N,n);
    S_0=100;
    Spath1 = S_0*[ones(N,1),cumprod(exp((r-0.5*sig^2)*dt + ...
                                       sig*sqrt(dt)*ranvec),2)];
    arithave1 = mean(Spath1,2);
    Parith1 = exp(-r*T)*max(arithave1-K,0); % payoffs

    %find f(x+h)  in the terms (f(x+h) - f(x))/h
    S_0 = S_0+h;
    Spath2 = S_0*[ones(N,1),cumprod(exp((r-0.5*sig^2)*dt + ...
                                   sig*sqrt(dt)*ranvec),2)];
     arithave2 = mean(Spath2,2);
    Parith2 = exp(-r*T)*max(arithave2-K,0); % payoffs

    delta=(Parith2 - Parith1)/h;

    D = mean(delta);
    stderr = std(delta)/sqrt(N);
    result3=[result3; K D stderr];
    K=K+10
end

disp ('     K        Delta(LH)  stderr');
result3

%Vega by All Methods}
%====================

% likelihood method for vega
clear
r=0.05; % annual interest
sig=0.3; % volatility
K=35; % stike price
S_0=100; % initial stock price
T=4/12; %maturity in 4 months, which is 4/12 of the year
n= 88;  % there are approx. 88 trading days in 4 months
dt=T/n; % time step
% Simulate N sample paths of the stock process

N=10^5;
K=20
result1=[];
for j=1:14
 %Vega=nan(N,1);
   for i=1:N
       Z=randn(1,n);
       path=(r-sig^2/2)*dt+sig*sqrt(dt)*Z;
       path=cumprod([S_0,exp(path)]);
       A=mean(path); % average path/stock price
       kk=((Z.^2-1)/sig-sqrt(dt)*Z);
       Vega(i)=exp(-r*T)*max(0,A-K).*sum(kk);
   end
   V=mean(Vega);
   stderr = std(Vega)/sqrt(N);
   result1=[result1; K V stderr];
   K=K+10
end

disp ('    K         Vega(LH)  stderr');
result1
%-----------------------
% Vega by Pathwise method
% Book10.pdf page 541 .... book code
clear
r=0.05; %annual interest
sig=0.3; %volatility
S_0=100; %Initial stock price
T=4/12; %Maturity in 4 months
n=88;  % trading days
dt=T/n;  % time steps

result2=[];
K=20  % strike price
for j=1:14
   % Simulate N sample paths of the stock process
   N=10^5; Vega = nan(N,1);
   for i=1:N
      path = (r - sig^2/2)*dt + sig * sqrt(dt) * randn(1,n);
      path = cumprod([S_0, exp(path)]);
      A = mean(path);
      %Delta(i) = exp(-r*T) * (A>=K) * A/S_0;
      Vega(i) = exp(-r*T) * (A>=K) * sum ( path.*(log(path/S_0) - (r + .5*sig^2)*[0:dt:T]))/(n+1)/sig;
   end
   V = mean(Vega);
   stderr = std(Vega)/sqrt(N);
   result2=[result2; K  V  stderr];
   K=K+10
end

disp ('       K        V(PW)    stderr');
result2
%-----------------
% Vega by Finite Differences Method
clear
r=.05; % annual interest
sig=0.3; % volatility
K=35; % stike price
S_0=100; % initial stock price
T=4/12; %maturity in 4 months, which is 4/12 of the year
n= 88;  % there are approx. 88 trading days in 4 months
dt=T/n; % time step
% Simulate N sample paths of the stock process

h=0.003;  % assume 1% of sig, to be used in formula (f(x+h) - f(x))/h
N=10000;
K=20
result3=[];
for K=20:10:150
    %find f(x)  in the terms (f(x+h) - f(x))/h
    ranvec= randn(N,n);
    sig=0.3;
    Spath1 = S_0*[ones(N,1),cumprod(exp((r-0.5*sig^2)*dt + ...
                                       sig*sqrt(dt)*ranvec),2)];
    arithave1 = mean(Spath1,2);
    Parith1 = exp(-r*T)*max(arithave1-K,0); % payoffs

    %find f(x+h)  in the terms (f(x+h) - f(x))/h
    sig=sig+h;
    Spath2 = S_0*[ones(N,1),cumprod(exp((r-0.5*sig^2)*dt + ...
                                   sig*sqrt(dt)*ranvec),2)];
     arithave2 = mean(Spath2,2);
    Parith2 = exp(-r*T)*max(arithave2-K,0); % payoffs
    Vega=(Parith2 - Parith1)/h;

    V = mean(Vega);
    stderr = std(Vega)/sqrt(N);
    result3=[result3; K V stderr];
    K=K+10
end

disp ('     K        Vega(LH)  stderr');
result3
%------------------------------