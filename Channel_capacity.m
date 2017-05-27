clear Q P M i j x Q1 Pa Pa1 Q2 ;
% mu=2;
% lambda1=1;
% m=1;
% Ent=zeros(length(lambda1),1)
% for k=1:length(lambda1)
%     lambda=lambda1(k)

mu1=0.5;
lambda=0.5;
m=1;
Ent=zeros(length(mu1),1)
for k=1:length(mu1)
    mu=mu1(k)
    
%Derive Equation 1    
% it is basically distribution of del t -del tau . Where at point 100
% represent del t -del tau =0, range -10 to 10

i=1
for j=-10:0.1:10
Q(i)=integral(@(x)pdf('InverseGaussian',x,mu,lambda).*pdf('InverseGaussian',x+j,mu,lambda),0,1000);
i=i+1;
end


Q1=Q;

% Normal distributed del tao
X=[0:0.1:5];
Pa1 = normpdf(X,m,0.33);


% Ent=0;
% for i=1:50
% Ent=Ent+0.1*sum(Pa1(i).*(-(Q(Q>0).*(log2(Q(Q>0))))));

% end
% Probality matrix p(X,Y)

% conditional matix for -5 to 5 p(Y|X)
matrix_Y_X=zeros(51,101);
for i=1:51
    matrix_Y_X(i,:)=fliplr(Q(49+i:149+i));
end
matrix_Y_X([1:10:end],[1:10:end])

% calculte h(Y|X)
% reference Mackey p 159
Ent1=0;
for i=1:50
Ent1=Ent1-0.1*Pa1(i)*0.1*sum(Q(50+i:150+i).*log2(Q(50+i:150+i)));
end
Ent1

% claculate h(Y)
for i=1:150
Q2=Q1(i:i+50);
M(i)= 0.1* sum(fliplr(Q2).*Pa1);
end

% plot(M)
Ent2=0.1*sum(-(M(M>0).*(log2(M(M>0)))))

Ent(k)=Ent2-Ent1
end
plot(mu1,Ent)
% sum(-(M(M>0).*(log2(M(M>0)))))-sum(-(Q(Q>0).*(log2(Q(Q>0)))))
