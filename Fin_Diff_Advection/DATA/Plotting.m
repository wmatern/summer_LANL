%Plots the values from calculated data
for i = [1,5,10]
    type = 'FEM';
    S = load([type,'_',num2str(i),'.dat'],'-mat');
    phi = S.phi(3:end-2,3:end-2);
    phi_analyt = S.phi_analyt;
    clear S;
    h = figure(i); hold on, plot(0:length(phi)-1,phi(:,34),'o'), plot(0:length(phi)-1,phi_analyt(:,34),'r')
    %legend([type,'-SUPG'],'Analytical'), title([type,'-SUPG ','at t = ',num2str(i)]), 
    axis([0,length(phi)-1,-.5,1.5]); set(gca,'FontSize',14);
    hold off
    saveas(h,[type,'_',num2str(i),'.png']);
end
close 1 5 10
for i = [1,5,10]
    type = 'Lax-Wendroff_noTVD';
    S = load([type,'_',num2str(i),'.dat'],'-mat');
    phi = S.phi(3:end-2,3:end-2);
    phi_analyt = S.phi_analyt;
    clear S;
    h = figure(i); hold on, plot(0:length(phi)-1,phi(:,34),'o'), plot(0:length(phi)-1,phi_analyt(:,34),'r')
    %legend('Lax-Wendroff-noTVD','Analytical'), title(['Lax-Wendroff-noTVD ','at t = ',num2str(i)]), 
    axis([0,length(phi)-1,-.5,1.5]); set(gca,'FontSize',14);
    hold off
    saveas(h,[type,'_',num2str(i),'.png']);
end
close 1 5 10
for i = [1,5,10]
    type = 'Lax-Wendroff_superbee';
    S = load([type,'_',num2str(i),'.dat'],'-mat');
    phi = S.phi(3:end-2,3:end-2);
    phi_analyt = S.phi_analyt;
    clear S;
    h = figure(i); hold on, plot(0:length(phi)-1,phi(:,34),'o'), plot(0:length(phi)-1,phi_analyt(:,34),'r')
    %legend('Lax-Wendroff-superbee','Analytical'), title(['Lax-Wendroff-superbee ','at t = ',num2str(i)]), 
    axis([0,length(phi)-1,-.5,1.5]); set(gca,'FontSize',14);
    hold off
    saveas(h,[type,'_',num2str(i),'.png']);
end
close 1 5 10
%%
for i = [1,5,10]
    type = 'MPDATA';
    S = load([type,num2str(i),'.dat'],'-mat');
    phi = S.phi(3:end-2,3:end-2);
    phi_analyt = S.phi_analyt;
    clear S;
    h = figure(i); hold on, plot(0:length(phi)-1,phi(:,34),'o'), plot(0:length(phi)-1,phi_analyt(:,34),'r')
    %legend(type,'Analytical'), title([type,' at t = ',num2str(i)]), 
    axis([0,length(phi)-1,-.5,1.5]); set(gca,'FontSize',14);
    hold off
    saveas(h,[type,'_',num2str(i),'.png']);
end
close 1 5 10
%%
for i = [1,5,10]
    type = 'MUSCL';
    S = load([type,'_',num2str(i),'.dat'],'-mat');
    phi = S.phi(3:end-2,3:end-2);
    phi_analyt = S.phi_analyt;
    clear S;
    h = figure(i); hold on, plot(0:length(phi)-1,phi(:,34),'o'), plot(0:length(phi)-1,phi_analyt(:,34),'r')
    %legend(type,'Analytical'), title([type,' at t = ',num2str(i)]), 
    axis([0,length(phi)-1,-.5,1.5]); set(gca,'FontSize',14);
    hold off
    saveas(h,[type,'_',num2str(i),'.png']);
end
close 1 5 10