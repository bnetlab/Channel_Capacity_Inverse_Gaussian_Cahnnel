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
% %del t -del tau=0 at index 101
% %  range -10 to 10
 
% % Source distribution : Normal distributed del tao
X=[0:0.1:5];
Pa1 = normpdf(X,m,0.33);
Pa=zeros(1,50);
Pa2=[Pa Pa1] ;

%path T->M->R
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

% reference Mackey p 159
EntC1=0;
for i=1:100
EntC1=EntC1-0.1*Pa2(i)*0.1*sum(Q10(i:100+i).*log2(Q10(i:100+i)));
end
EntC1

% claculate h(Y)
for i=1:101
Q12=Q10(i:i+100);
M10(i)= 0.1* sum(fliplr(Q12).*Pa2);
end

% plot(M)
EntC2=0.1*sum(-(M10(M10>0).*(log2(M10(M10>0)))));

EntC(indexM)=EntC2-EntC1

end

plot(mu1,EntC)
hold on





