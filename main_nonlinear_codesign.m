% Clear WorkSpace 
opengl hardware % Feel free to comment this line if there is no aliasing issue

%% Create a codesign object and perform two PIs
sys = CDSys;
[J1_save, sys_mpi] = sys.modifiedPI();
[J2_save, sys_spi] = sys.standardPI();

%% Print out a simulation summary
disp('======================================')
disp(['The initial cost is ', num2str(J1_save(1))])
disp(['The final   cost is ', num2str(J1_save(end))])
disp(['The not codesigned cost is ', num2str(J2_save(end))])
disp(['The cost is reduced by ', ...
      num2str((J2_save(end)-J1_save(end))/J2_save(end)*100), ...
     '%, compared with no codesign'])
%% Print a summary of the parameters
sys_mpi.parameterReport();
 
%% Create the figure for cost vs. iterations
figure(1)
iter = 10;
plot((1:iter),log10(J1_save(1:iter)),'r-o',(1:iter), ...
	    log10(J2_save(1:iter)), 'b-^','Linewidth',2)
legend('Modified Policy Iteration', 'Convential Policy Iteration')
title('Convential v Modified policy iteration', 'Fontsize', 14)
% legend('J_{compac} (Plant redesign under system equivalence)',
%' J_a (Standard H2/H_{\infty} mixed control design)')
xlabel('Number of iterations', 'Fontsize', 14)
ylabel('Log_{10}(J)', 'Fontsize', 14)

%% Simulated the two systems obtained from different approaches
[t0,y0]=sys_spi.simulation();
[t,y]=sys_mpi.simulation();

%% Plot the system trajectories
figure(2)
tf = 3;
subplot(121)
plot(t0,y0(:,1)+1,'b-.',t,y(:,1)+1,'r-','Linewidth',2);
legend('Without co-design', 'With co-design')
xlabel('Time (sec)', 'Fontsize', 12)
ylabel('Load position xL', 'Fontsize', 12)
axis([0 tf 0 1.4])
subplot(122)
for i = 1:length(y0)
 uwc0(i) = -1/2 * sys_spi.R\ sys_spi.B'* sys_spi.dPhi(y0(i,1:4))'*sys_spi.w;
end
for i = 1:length(y)
 uwc(i) = -1/2 * sys_mpi.R\ sys_mpi.B'* sys_mpi.dPhi(y(i,1:4))'*sys_mpi.w;
end
plot(t0,uwc0,'b-.',t,uwc,'r-','Linewidth',2)
%+1,'b-.',t,y(:,1)+1,'r-','Linewidth',2)
legend('Without co-design', 'With co-design')%,'x_2','x_3','x_4')
xlabel('Time (sec)', 'Fontsize', 12)
ylabel('Control input u', 'Fontsize', 12)

%% Get the value function
[xx,yy,zz0] = sys_spi.getValueFunctionSurface();
[~,~,zz]  = sys_mpi.getValueFunctionSurface();

%% Plot the value function
figure(3)
surf(xx,yy,zz0)
hold on
surf(xx,yy,zz)
xlabel('x1:=xL-yd','Fontsize',14)
ylabel('x3:=xB','Fontsize',14)
zlabel('V(x_1,0,x_3,0)','Fontsize',14)
axis([-1.5 1.5 -2 2 0 1500])
hold off
view(gca,[-25.5 30]);

% Create textarrow
annotation(gcf,'textarrow',[0.194642857142857 0.180357142857143],...
	[0.457142857142857 0.576190476190476],'String',{'With','co-design'},...
	'FontSize',12,...
	'FontName','Helvetica');

% Create textarrow
annotation(gcf,'textarrow',[0.226785714285714 0.193164933135215],...
	[0.769047619047619 0.723132969034608],'String',{'Without co-design'},...
	'FontSize',12,...
	'FontName','Helvetica');

title('Comparison of the value functions','Fontsize',14)

