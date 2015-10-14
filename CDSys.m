classdef CDSys
    properties
        mL		= 1.5; % Initial value for mL
        mB		= 20;
        dL		= 15;  % Initial value for dL
        kB	    = 15; %7.5;
        kBn		= 1.5;
        dB		= 0.5;      
    end
    
    properties (SetAccess = private)
        mLmin   = 1;
        mL0;
        mLmax   = 3;
        
        mBmin   = 15;
        mB0;
        mBmax   = 25;
        
        dLmin   = 10;
        dL0;
        dLmax   = 20;
        
        kBmin    = 5; %5;
        kB0;
        kBmax    = 20; %10;
        
        dBmin   = 0.1;
        dB0;
        dBmax   =   1;
        
        Q      =    diag([1000, 1, 1, 1]);
        R      =   1;           % weighting matrix R
        x0     =   [-1 0 0 0]';  % Initial condition
        k      =   [-2 0 0 0];  % Initial stabilziing control policy
        tf     =   3;           % Simulation time duration
        
        w;
    end
    
    properties (Dependent)
        A;
        B;
        p;
        w0;
        ll;
    end
    
    methods
        function obj = CDSys()
            %Initialization
            obj = reset_w(obj);
            obj.mL0 = obj.mL;
            obj.mB0 = obj.mB;
            obj.dL0 = obj.dL;
            obj.kB0 = obj.kB;
            obj.dB0 = obj.dB;
           
        end
        
        function A = get.A(obj)
            A = [0    1             0       0;
                0 -obj.dL/obj.mL-obj.dL/obj.mB  obj.kB/obj.mB obj.dB/obj.mB;
                0     0             0       1;
                0 obj.dL/obj.mB -obj.kB/obj.mB  -obj.dB/obj.mB];
        end
        
        function B = get.B(obj)
            B = [0;
                1/obj.mL+1/obj.mB;
                0
                -1/obj.mB];
        end
        
        function p = get.p(obj)
            p = [obj.mL; obj.mB; obj.dL; obj.kB; obj.dB];
        end
        
        function ll = get.ll(obj)
            ll = length(obj.dPhi([0,0,0,0]));
        end
        
        function w0 = get.w0(obj)
            w0   =   [0 0 0 0 2*obj.R/obj.B(2)*(-obj.k(1)) zeros(1,obj.ll-5)]';
        end
        
        function J = getCost(obj)
            ycost = obj.simulation4calculation();
            J =  ycost(end,end);
        end
        
        function obj = reset_w(obj)
            obj.w = obj.w0;
        end
        
        function dX = LP_model_calc(obj, t, X)
            x       = X(1:4);
            x = x(:);
            ul      =  obj.k * x; % linear portion of the controller
            unl     =  -1/2 * obj.R \ obj.B' *obj.dPhi(x)' * ...
                [zeros(10,1); obj.w(11:end)] ;
            dx      =  obj.A * x + [0;obj.kB/obj.mB*x(3)^3;0; ...
                -obj.kB/obj.mB*x(3)^3]+ obj.B * (ul + unl);
            dH      =  x * x';
            dHv     =  dH(:);
            df      =  x * obj.R * unl;
            dJ      =  x'*obj.Q*x + (ul+unl) * obj.R *(ul+unl);
            dX      = [dx ; dHv ; df ; dJ];
        end
        
        function y = dPhi(~,x)
            
            x1 = x(1);
            x2 = x(2);
            x3 = x(3);
            x4 = x(4);
            
            y  = [ 2*x1      0          0        0; %x1^2   #1
                0         2*x2       0        0; %x2^2   #2
                0         0        2*x3       0; %x3^2   #3
                0         0          0      2*x4;%x4^2   #4
                x2        x1         0        0; %x1*x2  #5
                x3        0          x1       0; %x1*x3  #6
                x4        0          0        x1;%x1*x4  #7
                0         x3        x2        0; %x2*x3  #8
                0         x4         0        x2;%x2*x4  #9
                0         0         x4       x3; %x3*x4  #10
                4*x1^3       0          0        0; %x1^4   #11
                0       4*x2^3       0        0; %x2^4   #12
                0         0       4*x3^3      0; %x3^4   #13
                0         0          0     4*x4^3; %x4^4       #14
                2*x1*x2^2  2*x2*x1^2    0        0; %x1^2*x2^2    #15
                0       2*x2*x3^2  2*x3*x2^2  0; %x2^2*x3^2    #16
                0         0        2*x3*x4^2  2*x4*x3^2; %x3^2*x4^2  #17
                2*x1*x4^2    0          0        2*x4*x1^2; %x1^2*x4^2  #18
                2*x1*x3^2    0        2*x3*x1^2  0; %x1^2*x3^2          #19
                0       2*x2*x4^2    0        2*x4*x2^2; %x1^2*x3^2  #20
                ];
        end
        
        function obj = oneStepGradientPolicyImprovement(obj)
            X = [];
            Y = [];
            for x1 = -1.5:.5:1.5
                for x2 = -1.5:.5:1.5
                    for x3 = -1.5:.5:1.5
                        for x4 = -2:0.5:2
                            dP = obj.dPhi([x1 x2 x3 x4]);
                            
                            u = - 1/2 * inv(obj.R) * obj.B' * ...
                                obj.dPhi([x1 x2 x3 x4])' * obj.w;
                            
                            X = [X;(dP*(obj.A*[x1;x2;x3;x4] + ...
                                [0;obj.kBn/obj.mB*x3^3;0; ...
                                -obj.kBn/obj.mB*x3^3] + obj.B*u))'];
                            
                            Y = [Y; [x1;x2;x3;x4]'*obj.Q*[x1;x2;x3;x4]+ ...
                                u'*obj.R*u];
                            
                        end
                    end
                end
            end
            obj.w = X\(-Y); % update the weights
            obj.k = -1/2*inv(obj.R) * ...
                [obj.B(2)*obj.w(5) + obj.B(4)*obj.w(7) ...
                2*obj.B(2)*obj.w(2)+obj.B(4)*obj.w(9) ...
                obj.B(2)*obj.w(8)+obj.B(4)*obj.w(10) ...
                2*obj.B(4)*obj.w(4)+obj.B(2)*obj.w(9)];
        end
        
        function obj = oneStepParametersImprovement(obj)
            %    preparing the equality constraints for QP solvers
            %          k(1) k(2) k(3) k(4) mL       mB       dL   kB  dB
            Aeq = [ 1  0  0  0  -obj.k(1)/obj.mL          0      0   0   0 %k1/mL
                0  1  0  0  -(obj.k(2)-obj.dL)/obj.mL     0     -1   0   0 %k2/mL-dL/mL
                0  0  1  0  -obj.k(3)/obj.mL          0      0   0   0 %k3/mL
                0  0  0  1  -obj.k(4)/obj.mL          0      0   0   0 %k4/mL
                1  0  0  0      0         -obj.k(1)/obj.mB   0   0   0 %-k1/mB
                0 -1  0  0      0    -(obj.dL-obj.k(2))/obj.mB   1   0   0 %dL/mB-k2/mB
                0  0 -1  0      0     (obj.k(3)+obj.kB)/obj.mB   0  -1   0 %-k3/mB-kB/mB
                0  0  0 -1      0     (obj.k(4)+obj.dB)/obj.mB   0   0  -1 ];%-k4/mB-dB/mB
            Beq = zeros(8,1);
            
            
            
            [Aeq, Beq] = obj.mat_rdc(Aeq, Beq); % extract the independent rows from Aeq
            
            %preparing inequality constraints for QP solvers
            %        k1    k2  k3  k4    mL     mB    dL    kB    dB
            LB  = [-Inf -Inf -Inf -Inf  obj.mLmin obj.mBmin obj.dLmin obj.kBmin obj.dBmin]';
            UB  = [ Inf  Inf  Inf  Inf  obj.mLmax obj.mBmax obj.dLmax obj.kBmax obj.dBmax]';
            PK=[obj.k(:);obj.p];
            
            
            % 			[~,yc] = ode45(@(t,X) obj.LP_model_calc(t,X), ...
            % 				[0, 100],[obj.x0;zeros(4*4+4+1,1)]);
            % 			H       =  obj.R * reshape(yc(end, 4+1:4+16)',4,4);
            % 			f       =  yc(end, 4+16+1:4+16+4)';
            
            [yc, H, f] = obj.simulation4calculation();
            
            options.LargeScale = 'off';  % setting opitons for the Q_Programming
            options.Display = 'off';     % get rid of those super annoying text
            
            PK = quadprog(blkdiag(H,zeros(5)),[f;zeros(5,1)],[],[], ...
                Aeq, Beq, LB, UB, PK, options);
            
            obj.k =  PK(1:4)'; %update the linear part of the controller
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Update the co-design parameters
            
            obj.mL = PK(5);
            obj.mB = PK(6);
            obj.dL = PK(7);
            obj.kB = PK(8);
            obj.dB = PK(9);
        end
        
        function [J1_save, obj_mpi] = modifiedPI(obj)
            tic()
            disp('Working on the modified PI...')
            J1_save = obj.getCost();
            for i=1:10
                disp([num2str(i), '-th iteration'])
                obj = oneStepGradientPolicyImprovement(obj);
                obj = oneStepParametersImprovement(obj);
                J1_save = [J1_save, obj.getCost()];
            end
            modified_PI_cpu = toc()
            obj_mpi = obj;
        end
        
        function [J2_save, obj_spi] = standardPI(obj)
            tic()
            disp('Working on the standard PI...')
            J2_save = obj.getCost();
            for i=1:10
                disp([num2str(i), '-th iteration'])
                obj = oneStepGradientPolicyImprovement(obj);
                J2_save = [J2_save, obj.getCost()];
            end
            stand_PI_cpu = toc()
            obj_spi = obj;
        end
        
        function [xx,yy,zz] = getValueFunctionSurface(obj)
            figure(3)
            xx = -1.5:0.15:1.5;
            yy = -1.5:0.15:1.5;
            zz = zeros(length(yy), length(xx)); 
            
            for i = 1:length(xx)
                for j = 1:length(yy)
                    zz(j,i) = obj.Phi([xx(i),0,yy(j),0])*obj.w;
                end
            end
            
            [xx,yy]=meshgrid(xx,yy);
        end
        
        function [tout,yout] = simulation(obj)
            [tout,yout] = ode45(@(t,X) obj.LP_model(t,X),[0,obj.tf],obj.x0);
        end
        
        function [yc, H,f] = simulation4calculation(obj)
            [~,yc] = ode45(@(t,X) obj.LP_model_calc(t,X), ...
                [0, 100],[obj.x0;zeros(4*4+4+1,1)]);
            H       =  obj.R * reshape(yc(end, 4+1:4+16)',4,4);
            f       =  yc(end, 4+16+1:4+16+4)';
        end
        
        function dx = LP_model(obj,t,x)
            ul       =  obj.k * (x(1:4)); % linear portion of the controller
            unl      =  -1/2*inv(obj.R) * obj.B' *obj.dPhi(x)' * [zeros(10,1); obj.w(11:end)] ;
            dx       =  obj.A * x + [0;obj.kB/obj.mB*x(3)^3;0;-obj.kB/obj.mB*x(3)^3]+ obj.B * (ul + unl);
        end
        
        function [An,Bn] = mat_rdc(~,A,B)
            n=length(B);
            i=0;
            An=[];
            Bn=[];
            i=0;
            while i < n
                i = i+1;
                if rank([An; A(i,:)])>rank(An)
                    An = [An; A(i,:)];
                    Bn = [Bn; B(i,:)];
                end
            end
        end
        
        function parameterReport(obj)
            disp('======================================')
            fprintf('Variable Min  Max   Initial  Optimized\n' )
            fprintf('mL       %4.2f  %4.2f  %4.2f     %4.2f\n', ...
                obj.mLmin, obj.mLmax, obj.mL0, obj.mL )
            fprintf('mB      %4.2f %4.2f %4.2f    %4.2f\n', ...
                obj.mBmin, obj.mBmax, obj.mB0, obj.mB )
            fprintf('dL      %4.2f %4.2f %4.2f    %4.2f\n', ...
                obj.dLmin, obj.dLmax, obj.dL0, obj.dL )
            fprintf('kB       %4.2f %4.2f %4.2f    %4.2f\n', ...
                obj.kBmin, obj.kBmax, obj.kB0, obj.kB )
            fprintf('dB       %4.2f  %4.2f  %4.2f     %4.2f\n', ...
                obj.dBmin, obj.dBmax, obj.dB0, obj.dB )
            disp('======================================')
        end
        
        function y = Phi(~,x)
            
            x1 = x(1);
            x2 = x(2);
            x3 = x(3);
            x4 = x(4);
            
            y  = [x1^2
                x2^2
                x3^2
                x4^2
                x1*x2
                x1*x3
                x1*x4
                x2*x3
                x2*x4
                x3*x4
                x1^4
                x2^4
                x3^4
                x4^4
                x1^2*x2^2
                x2^2*x3^2
                x3^2*x4^2
                x1^2*x4^2
                x1^2*x3^2
                x1^2*x3^2
                ]';
        end
    end
end