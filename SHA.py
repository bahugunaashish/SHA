load a4.txt
x3 = a4(1,1); 
y3 = a4(1,2);

  for col=1:17
      load f1D.txt
     x = f1D(col,1:2); 
     y = f1D(col,3:4);
     M_max=f1D(col,5);
     b=f1D(col,6);
     a=f1D(col,7);
     L= ((y(1)-y(2))^2+(x(1)-x(2))^2)^0.5;
     M = (y(2)-y(1))/(x(2)-x(1));
     X4= ((y3-y(1))+1/M*x3+M*x(1))/(M+1/M);
     Y4= -(X4-x3)/M+y3;
 
check1=((Y4-y(1))^2+(X4-x(1))^2)^0.5;
check2=((Y4-y(2))^2+(X4-x(2))^2)^0.5;
check=max(check1,check2);
if (check<=L)
R_MIN = ((X4-x3)^2+(Y4-y3)^2)^0.5;
R1 = ((x(1)-x3)^2+(y(1)-y3)^2)^0.5;
R2 = ((x(2)-x3)^2+(y(2)-y3)^2)^0.5;
R_MAX = max(R1,R2);
end
if(L<check)
R1 = ((x(1)-x3)^2+(y(1)-y3)^2)^0.5;
R2 = ((x(2)-x3)^2+(y(2)-y3)^2)^0.5;
R_MIN= min(R1,R2);
R_MAX = max(R1,R2);
end
T = 10;
D = (R_MAX-R_MIN)/T;
Z1 = (R1^2-R_MIN^2)^0.5;
Z2 =(R2^2-R_MIN^2)^0.5;
perpnd=((X4-x3)^2+(Y4-y3)^2)^0.5;
    extra=(R_MIN^2-perpnd^2)^0.5;
    
if(L<check)
 
    i=1;
    N=1;
     while ((R_MIN+N*D)<=R_MAX)
           DS(i)=((R_MIN+i*D)^2-perpnd^2)^0.5;
           N=N+1;
           DiS(i)=DS(i);
           i=i+1;
    end
    for i=1:i-1
       % DiS(i);
       % fprintf('\nDiS= %.2f\n',DiS(i));
        
    end
  
s=1;
counter=0;
while (s<i+1)
if (s==1)
    DiS(s)=DS(s)-extra;
else
    DiS(s)=DiS(s)-DS(s-1);
end

%fprintf('distance2(%d)=%f\n\n\n',s,DiS(s));
    
s=s+1;
end
for i=1:T
     dt(i)=DiS(i)/L;
     %fprintf('P(r)(%d)=%f\n\n',i,dt(i));
end
r=linspace(R_MIN,R_MAX,T+1);
 for i=1:T
     if (i<=T)
     R_avg(i)=(r(i)+r(i+1))/2;
     
     end
 end
     figure (1);
 bar(R_avg,dt)
end

if (check<=L)
i=1;
N=1;
while ((R_MIN+N*D)<=R1)
    DS(i)=((R_MIN+i*D)^2-R_MIN^2)^0.5;
    N=N+1;
    DiS(i)=DS(i);
    i=i+1;
end
for i=1:i-1
    %DiS(i);
    %fprintf('\nDiS = %.2f\n',DiS(i));
end

s=1;
counter=0;
while (s<i+1)
if (s==1)
    DiS(s)=DS(s);
else
    DiS(s)=DiS(s)-DS(s-1);
end


counter=counter+DiS(s);    
s=s+1;
end
X=s;
    counter;
remain1=Z1-counter;

i=1;
n=1;
while ((R_MIN+n*D)<=R2)
    DS1(i)=((R_MIN+i*D)^2-R_MIN^2)^0.5;
    n=n+1;
    DiS1(i)=DS1(i);
    i=i+1;
end
for i=1:i-1
    DiS1(i);
    %fprintf('\nDis1= %.2f\n',DiS1(i));
end

s=1;
counter=0;
while (s<i+1)
if (s==1)
    DiS1(s)=DS1(s);
else
    DiS1(s)=DiS1(s)-DS1(s-1);
end

counter=counter+DiS1(s);    
s=s+1;
end
sit=s;
remain2=Z2-counter;

i=1;
n=1;
if(R1>R2)
while(i<=T)
    if(n==s)
            dt(i)=(DiS(n)+remain2)/L;
    end
    if(n<s) 
    dt(i)=(DiS(i)+DiS1(n))/L;
    end
    if(n>s)
     dt(i)=DiS(i)/L;
    end
   % fprintf('P(r)(%d)=%d\n\n',i,dt(i));
    i=i+1;
    n=n+1;
 
end
end
i=1;
n=1;
if (R2>R1)
    while(i<=T)
    if(n==X)
            dt(i)=(DiS1(n)+remain1)/L;
    end
    if(n<X) 
       dt(i)=(DiS1(i)+DiS(n))/L;
    end
    if(n>X)
      dt(i)=DiS1(i)/L;
    end
  
    i=i+1;
    n=n+1;
end
 
end

 r=linspace(R_MIN,R_MAX,T+1);
 for i=1:T
     if (i<=T)
     R_avg(i)=(r(i)+r(i+1))/2;
     
     end
    
 end
 figure (2);
 bar(R_avg,dt)
end 
 
 
 m0=4;
beta = 2.303*b;
interval=(M_max-m0)/T;
i=1;
minimum=m0;
while (i<=T+1)
    if (i==1)
    mx(i)=minimum;
    end
    if (i>1)
    mx(i)=minimum+interval;
    minimum=mx(i);
    end
    i=i+1;
end
s=1;
 
while(s<=T)
    pm(s)=(exp(-beta*(mx(s)-m0))-exp(-beta*(mx(s+1)-m0)))/(1-exp(-beta*(M_max-m0)));
    s=s+1;
end
h=1;
for h=1:T
     if (h<=T)
     M_avg(h)=(mx(h)+mx(h+1))/2;
     end
end
figure (5);
bar(M_avg,pm)

 
pha=[0.01 3];
SIGMA_LN=0.57;
cx=pha(2)/0.01;
for i=1:cx
    if(i==1)
PHA(i)=pha(i);
    end
    if(i>1)
PHA(i)=PHA(i-1)+0.01;
    end
end

yx=0;
mvx=1;  
i=1;
while (mvx<=T)
             rvx=1;
          while (rvx<=T)
              
                %tempo(i)=exp(3.374+0.3503*(M_avg(mvx)-6)-0.0698*(M_avg(mvx)-6)^2-log(R_avg(rvx))-0.00919*R_avg(rvx))/981;
                %LN_PHA(i)=log(tempo(i));
             
               LN_PHA(i) = -6.680+1.134*M_avg(mvx)-0.001*R_avg(rvx)-0.7098*log(R_avg(rvx));  
                %LN_PHA(i) = 6.74+0.859*M_avg(mvx)-1.80*log(R_avg(rvx)+25); 
              
                 rvx=rvx+1;
                 i=i+1;
          end
        mvx=mvx+1;
end
ux=i;
w=1;
count=1;
while(w<=cx)
            LQ=1;
            while (LQ<=ux-1)
                Z_value(count)=(log(PHA(w))-LN_PHA(LQ))/SIGMA_LN;
              LQ=LQ+1;
            count=count+1;
            end
        w=w+1;
  
end
i=1;
while(i<=(cx*T^2))
     pz_value(i)=1-normcdf(Z_value(i),0,1);
     i=i+1;
end

c=1;
z=1;
while(z<=cx)
    i=1;
      while (i<=10)
             j=1;
             while(j<=10)
                 k(c)=pm(i)*dt(j)*pz_value(c);
                   j=j+1;
                   c=c+1;
             end
             i=i+1;
      end
       z=z+1;
end
 
 

lemda=10^(a-b*m0);
k=lemda*k;
 

 
i=1;
c=1;
while(i<=cx)
    annual_exced=0;
    while(c<(i*T^2))
        annual_exced=annual_exced+k(c);
        c=c+1;
    end
    anexc(i)=annual_exced;
    i=i+1;
end

  if (col==1)
     anexc_a_3=anexc;
     anexc_col_a_3=anexc';
 
 semilogy(PHA,anexc_a_3,'r');
 end
 if (col==2)
     anexc_ccf=anexc;
     anexc_col_ccf=anexc';
 
 semilogy(PHA,anexc_ccf,'r');
 end
 if (col==3)
     anexc_cmf=anexc;
     anexc_col_cmf=anexc';
 
 semilogy(PHA,anexc_cmf,'r');
 end
 if (col==4)
     anexc_dauk=anexc;
     anexc_col_dauk=anexc';
 
 semilogy(PHA,anexc_dauk,'r');
 end
 if (col==5)
     anexc_dhub=anexc;
     anexc_col_dhub=anexc';
 
 semilogy(PHA,anexc_dhub,'r');
 end
 if (col==6)
     anexc_naga=anexc;
     anexc_col_naga=anexc';
 
 semilogy(PHA,anexc_naga,'r');
 end
 if (col==7)
     anexc_ebt_k=anexc;
     anexc_col_ebt_k=anexc';
 
 semilogy(PHA,anexc_ebt_k,'r');
 end
 if (col==8)
     anexc_kalad=anexc;
     anexc_col_kalad=anexc';
 
 semilogy(PHA,anexc_kalad,'r');
 end
 if (col==9)
     anexc_kopil=anexc;
     anexc_col_kopil=anexc';
 
 semilogy(PHA,anexc_kopil,'r');
 end
 if (col==10)
     anexc_lohi=anexc;
     anexc_col_lohi=anexc';
 
 semilogy(PHA,anexc_lohi,'r');
 end
 if (col==11)
     anexc_mbt_mct=anexc;
     anexc_col_mbt_mct=anexc';
 
 semilogy(PHA,anexc_mbt_mct,'r');
 end
 if (col==12)
     anexc_mish=anexc;
     anexc_col_mish=anexc';
 
 semilogy(PHA,anexc_mish,'r');
 end
 if (col==13)
     anexc_cb=anexc;
     anexc_col_cb=anexc';
 
 semilogy(PHA,anexc_cb,'r');
 end
 if (col==14)
     anexc_old=anexc;
     anexc_col_old=anexc';
 
 semilogy(PHA,anexc_old,'r');
 end
 if (col==15)
     anexc_sag=anexc;
     anexc_col_sag=anexc';
 
 semilogy(PHA,anexc_sag,'r');
 end
 if (col==16)
     anexc_sam=anexc;
     anexc_col_sam=anexc';
 
 semilogy(PHA,anexc_sam,'r');
 end
 if (col==17)
     anexc_sylh=anexc;
     anexc_col_sylh=anexc';
 
 semilogy(PHA,anexc_sylh,'r');
 end

    end

%figure(7); 
 %semilogy(PHA,anexc_a_3,'r');

%hold on 
%semilogy(PHA,anexc_ccf,'r');

 
 %semilogy(PHA,anexc_cmf,'r');

 
 %semilogy(PHA,anexc_dauk,'r');

 
 %semilogy(PHA,anexc_dhub,'r');

 
 %semilogy(PHA,anexc_naga,'r');

 
 %semilogy(PHA,anexc_ebt_k,'r');

 
 %semilogy(PHA,anexc_kalad,'r');

 
 %semilogy(PHA,anexc_kopil,'r');


 
 %semilogy(PHA,anexc_lohi,'r');

 
 %semilogy(PHA,anexc_mbt_mct,'r');

 
% semilogy(PHA,anexc_mish,'r');

 
% semilogy(PHA,anexc_cb,'r');

% semilogy(PHA,anexc_old,'r');


 
 %semilogy(PHA,anexc_sag,'r');


 
 %semilogy(PHA,anexc_sam,'r');
 
 
 
 %semilogy(PHA,anexc_sylh,'r');
 %hold off



