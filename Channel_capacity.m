clear Q P M i j x Q1 Pa Pa1 Q2 Q10 EntC Ent ;
clear all;
% parameter value and looping
mu=0.3;
lambda1=[0.1 0.2 0.3 0.4 0.5 1];
m=1;
Ent=zeros(length(lambda1),1)
for indexM=1:length(lambda1)
    lambda=lambda1(indexM)

% mu1=[0.1 0.15 0.2 0.3 0.4 0.5 1];
% lambda=0.5;
% m=1;
% Ent=zeros(length(mu1),1)
% for indexM=1:length(mu1)
%     mu=mu1(indexM)
    
% %Derive Equation 1    
% % it is basically distribution of del t -del tau . 
% % del t -del tau=0 at index 101
% %  range -10 to 10

i=1;
disp('1st scenario')
for j=-10:0.1:10
Q(i)=integral(@(x)pdf('InverseGaussian',x,mu,lambda).*pdf('InverseGaussian',x+j,mu,lambda),0,1000);
i=i+1;
end

Q1=Q;
Q13=Q;
 
% % Source distribution : Normal distributed del tao
X=[0:0.1:5];
Pa1 = normpdf(X,m,0.33);
Pa=zeros(1,50);
Pa2=[Pa Pa1] ;

% Probality matrix p(X,Y) of the channel

% conditional matix for -5 to 5 p(Y|X)
% matrix_Y_X=zeros(51,101);
% for i=1:101
%     matrix_Y_X(i,:)=fliplr(Q(i:100+i));
% end
% matrix_Y_X([1:10:end],[1:10:end])


% %Entropy calculation

% Ent=0;
% for i=1:50
% Ent=Ent+0.1*sum(Pa1(i).*(-(Q(Q>0).*(log2(Q(Q>0))))));
% end

% % calculte h(Y|X)
% % reference Mackey p 159
% Ent1=0;
% for i=1:100
% Ent1=Ent1-0.1*Pa2(i)*0.1*sum(Q(i:100+i).*log2(Q(i:100+i)));
% end
% Ent1
% 
% % claculate h(Y)
% for i=1:101
% Q2=Q1(i:i+100);
% M(i)= 0.1* sum(fliplr(Q2).*Pa2);
% end
% 
% % plot(M)
% Ent2=0.1*sum(-(M(M>0).*(log2(M(M>0)))))
% 
%  Ent(k)=Ent2-Ent1
% sum(-(M(M>0).*(log2(M(M>0)))))-sum(-(Q(Q>0).*(log2(Q(Q>0)))))

% % Expresion for indirect path T->M->R

i=1;
disp('2st scenario')
for j=-10:0.1:10
QY(i)=integral(@(x)pdf('InverseGaussian',x,mu/2,lambda/4).*pdf('InverseGaussian',x+j,mu/2,lambda/4),0,1000);
i=i+1;
end

X=[-10:.1:10];
ind=1;
clear Q10
for value = -10:0.1:0
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

Q23=Q10;

% reference Mackey p 159
% % EntC1=0;
% % for i=1:100
% % EntC1=EntC1-0.1*Pa2(i)*0.1*sum(Q10(i:100+i).*log2(Q10(i:100+i)));
% % end
% % EntC1
% % 
% % % claculate h(Y)
% % for i=1:101
% % Q12=Q10(i:i+100);
% % M10(i)= 0.1* sum(fliplr(Q12).*Pa2);
% % end
% % 
% % % plot(M)
% % EntC2=0.1*sum(-(M10(M10>0).*(log2(M10(M10>0)))));
% % 
% % EntC(k)=EntC2-EntC1
% % trapz(Q10)
% % end

% % plot(lambda1,EntC)
% % hold on
% % plot(lambda1,Ent)


% % 3rd scenario

disp('3rd scenario')
i=1
j=0
for T=-10:0.1:10
Q30(i)=integral(@(x)pdf('InverseGaussian',x-j,mu/2,lambda/4).*pdf('InverseGaussian',T-x,mu/2,lambda/4),j,1000);
i=i+1;
end

count=1;
for k=-10:0.1:10
    i=1;
    for T=-10:0.1:10
        Q31(i)=pdf('InverseGaussian',T-k,mu,lambda);
        i=i+1;
    end
Q32=Q30.*Q31;
Q33(count)=0.1*trapz(Q32);
count=count+1;
end
trapz(Q33)

% % 4th scenario
disp('4th scenario')
i=1;
for T=-10:0.1:10
Q40(i)=integral(@(x)pdf('InverseGaussian',x,mu/2,lambda/4).*pdf('InverseGaussian',T-x,mu/2,lambda/4),0,1000);
i=i+1;
end

count=1;
j=0;
for k=-10:0.1:10
    i=1;
    for T=-10:0.1:10
        Q41(i)=pdf('InverseGaussian',T+k-j,mu,lambda);
        i=i+1;
    end
Q42=Q40.*Q41;
Q43(count)=0.1*trapz(Q42);
count=count+1;
end
trapz(Q43)
% %normalization

Q50=Q13+Q23+Q33+Q43;
Q50=Q50/4;


% % Mutual information
Ent1=0;
for i=1:100
Ent1=Ent1-0.1*Pa2(i)*0.1*sum(Q50(i:100+i).*log2(Q50(i:100+i)));
end


% % claculate h(Y)
for i=1:101
Q2=Q50(i:i+100);
M(i)= 0.1* sum(fliplr(Q2).*Pa2);
end

% % plot(M)
Ent2=0.1*sum(-(M(M>0).*(log2(M(M>0)))));

Ent(indexM)=Ent2-Ent1
end
plot(lambda1,Ent)