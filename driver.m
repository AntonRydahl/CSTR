function driver(NS)
close all;
% reading C output

time_steps_per_sample = 60;
number_of_samples = 35;
total_steps = number_of_samples*time_steps_per_sample;
N = total_steps;
n=3;

fileID = fopen('X.txt','r');
formatSpec = '%f';
A = fscanf(fileID,formatSpec);
X = reshape(A,[n,(N+1),NS]);

fileID2 = fopen('T.txt','r');
T = fscanf(fileID2,formatSpec);

fileID2 = fopen('F.txt','r');
F = fscanf(fileID2,formatSpec);

u = @(t,F)(F(floor(t)+1));

% F = fscanf(fileID3,formatSpec);

% Plotting output

dim = [2 2];
figure('Renderer', 'painters', 'Position', [10 10 1500 1000])
subplot(dim(1),dim(2),1)
plot(T,X(1,:,1))
hold on
for i=2:NS
    plot(T,X(1,:,i))
end
hold off
xlabel('$t\:[min]$','interpreter','latex','fontsize',16)
ylabel('$C_A [\frac{mol}{L}]$','interpreter','latex','fontsize',16)

subplot(dim(1),dim(2),2)
plot(T,X(2,:,1))
hold on
for i=2:NS
    plot(T,X(2,:,i))
end
hold off
xlabel('$t\:[min]$','interpreter','latex','fontsize',16)
ylabel('$C_B [\frac{mol}{L}]$','interpreter','latex','fontsize',16)

subplot(dim(1),dim(2),3)
plot(T,X(3,:,1))
hold on
for i=2:NS
    plot(T,X(3,:,i))
end
hold off
xlabel('$t\:[min]$','interpreter','latex','fontsize',16)
ylabel('$T [K]$','interpreter','latex','fontsize',16)

subplot(dim(1),dim(2),4)
plot(T(1:end-1),u(T(1:end-1),F),'Color','#FF8C00','linewidth',2)
xlabel('$t\:[min]$','interpreter','latex','fontsize',16)
ylabel('$F\:[\frac{mL}{min}]$','interpreter','latex','fontsize',16)
ylim([0 1000])

solverName = 'implicit_explicit';

title_str="Stochastic Solution to the CSTR Model";

sgtitle(title_str,'fontsize',22)
set(gca,'LooseInset',get(gca,'TightInset'));
savefig(gcf,strcat([solverName,'.fig']))
saveas(gcf,strcat([solverName,'.png']))

end

