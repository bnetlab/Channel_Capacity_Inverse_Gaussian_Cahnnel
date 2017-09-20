clear all;
% % parameter value and looping

% V=2;
% sigma=0.5;
% d=20;
% mu= d./V
% lambda1=d.^2./sigma.^2
% 
% Ent=zeros(length(lambda1),1)
% for indexM=1:length(lambda1)
%     lambda=lambda1(indexM)

V=[2 3 1.5];
sigma=[0.2 0.3 0.5];
d=[10 20 30];
count1=1;

for ii=1:3
    for jj=1:3
        for kk=1:3

mu= d(kk)./V(ii)
lambda=d(kk).^2./sigma(jj).^2

% mu1=[0.33 0.40 0.50 0.667 1.00 2.00];
% lambda=1;
% m=1;
% Ent=zeros(length(mu1),1)
% for indexM=1:length(mu1)
%     mu=mu1(indexM)
    
% %Derive Equation 1    
% % it is basically distribution of del t -del tau . 
% %del t -del tau=0 at index 101
% %  range -10 to 10
 
% % Source distribution : Normal distributed del tao
% X=[0:0.1:5];
% Pa1 = normpdf(X,0.5,0.25);
% Pa=zeros(1,50);
% Pa2=[Pa Pa1] ;

% % uniform
X=[0:0.1:5];
Pa1=unifpdf(X,0,2);
Pa=zeros(1,50);
Pa2=[Pa Pa1];

% % %exp
% X=[0:0.1:5];
% Pa1=exppdf(X,1.359);
% Pa=zeros(1,50);
% Pa2=[Pa,Pa1];

%path T->M->R
i=1;
disp('2st scenario')
for j=-10:0.1:10
QY(i)=integral(@(x)pdf('InverseGaussian',x,mu/2,lambda/4).*pdf('InverseGaussian',x+j,mu/2,lambda/4),0,1000);
i=i+1;
end

% Q10=0.1*conv(QY,QY);
% Q10=Q10(101:301);
thes=find(QY<0.0005);
QY(thes)=0;

X=[-100:1:100];
ind=1;
clear Q10
for value = -100:1:0
    count=1;
    for i=1:201 %%checking for sum of time deviation in both path is del t
        for j=1:201
            if (X(i)+X(j)==value)
                Xa(count)=i;
                Ya(count)=j;
                count=count+1;
            end
        end
    end
    Q10(ind)=sum(QY(Xa).*QY(Ya));
    ind=ind+1;
    clear Xa
    clear Ya
end

size(Q10)
Q10= 0.1*[Q10 fliplr(Q10(1:100))];
trapz(Q10)

% plot(Q10)
% hold on

% reference Mackey p 159
EntC1=0;
for i=1:100
EntC1=EntC1-0.1*Pa2(i)*0.1*sum(Q10(i:100+i).*mylog2(Q10(i:100+i)));
end
EntC1

% claculate h(Y)
for i=1:101
Q12=Q10(i:i+100);
M10(i)= 0.1* sum(fliplr(Q12).*Pa2);
end

% plot(M)
EntC2=0.1*sum(-(M10(M10>0).*(mylog2(M10(M10>0)))));

EntC(count1)=EntC2-EntC1

% 
% plot(lambda1,EntC)
% hold on

X=[-10:0.1:10];
avg=0.1*trapz(X.*Q10);
X1=(X-avg).^2;
Ycon(count1)=0.1*trapz(X1.*Q10)

% Legend=cell(6,1)%  two positions 
% Legend{1}=' 0.16' ;
% Legend{2}=' 0.25';
% Legend{3}=' 0.44' ;
% Legend{4}=' 1';
% Legend{5}=' 4' ;
% Legend{6}='16';
% legend(Legend);

AV(count1)=V(ii);
Asigma(count1)=sigma(jj);
Ad(count1)=d(kk);

count1=count1+1

        end
    end
end



