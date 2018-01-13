function [ML, MR] = calcMatrix(pTask, pCalc)
%%
hbar  = pTask.hbar;
B_inf = pTask.B_inf;
n     = pCalc.n;
h     = pCalc.h;      
tau   = pCalc.tau;
%%
switch pCalc.name
    case 'FFDS'
        theta = pCalc.theta;
        b = zeros(1, n); c = b;
        jj=1:n;
        b(jj) = hbar^2/(2*h)*B(pTask, h*(jj-3/2));
        c(jj) = h*(V(pTask, h*(jj-3/2))-1i*2*hbar/tau*rho(pTask, h*(jj-3/2)));
        a = -b+theta*c;
        %%
        Va = zeros(1, n-1); VaR = Va; % subdiagonal
        Vd = zeros(1, n  ); VdR = Vd; % diagonal
        Vc = Va;            VcR = Vc; % superdiagonal
        if pCalc.pLBoundary.bUseDTBC % left boundary
            R = pCalc.pLBoundary.R;
            c_0 = pCalc.pLBoundary.c_0;
            Vd(1)  = b(1)+(1/2-pCalc.pLBoundary.theta)*c(1)      -hbar^2*B_inf*c_0*R(1);
            VdR(1) = b(1)+(1/2-pCalc.pLBoundary.theta)*conj(c(1))-hbar^2*B_inf*c_0*R(2);
            Vc(1)  =-b(1)+     pCalc.pLBoundary.theta *c(1);
            VcR(1) = conj(Vc(1));
        else
            Vd(1) = 1; VdR(1) = 0; Vc(1) = 0; VcR(1) = 0;
        end;
        % inner nodes
        jj=2:(n-1);
        Va(jj-1)  = a(jj);
        Vd(jj)    = b(jj)+b(jj+1)+(1/2-theta)*(c(jj)+c(jj+1));
        Vc(jj)    = a(jj+1);
        VaR(jj-1) = conj(Va(jj-1));
        VdR(jj)   = conj(Vd(jj));
        VcR(jj)   = conj(Vc(jj));
        if pCalc.pRBoundary.bUseDTBC % right boundary
            R = pCalc.pRBoundary.R;
            c_0 = pCalc.pRBoundary.c_0;
            Vd(n)    = b(n)+(1/2-pCalc.pRBoundary.theta)*c(n)      -hbar^2*B_inf*c_0*R(1);
            VdR(n)   = b(n)+(1/2-pCalc.pRBoundary.theta)*conj(c(n))-hbar^2*B_inf*c_0*R(2);
            Va(n-1)  =-b(n)+     pCalc.pRBoundary.theta *c(n);
            VaR(n-1) = conj(Va(n-1));
        else
            Vd(n) = 1; VdR(n) = 0; Va(n-1) = 0; VaR(n-1) = 0;
        end;
        ML = gallery('tridiag', Va,  Vd,  Vc ); % full(gallery('tridiag', Va,  Vd,  Vc ));
        MR = gallery('tridiag', VaR, VdR, VcR); % full(gallery('tridiag', VaR, VdR, VcR));
%         keyboard
%         b1 = zeros(1, n-1); c1 = b1; % (n-1) - number of blocks in ML 
%         jj=1:(n-1);
%         b1(jj) = hbar^2/(2*h)*B(pTask, h*(jj-1/2));
%         c1(jj) = h*(V(pTask, h*(jj-1/2))-1i*2*hbar/tau*rho(pTask, h*(jj-1/2)));
%         ML1 = zeros(n);
%         for j=1:(n-1)
%             v = j:(j+1);
%             ML1(v, v) = ML1(v, v) + c1(j)*[1-2*theta 2*theta; 2*theta 1-2*theta]/2 + b1(j)*2*[1/2 -1/2; -1/2 1/2];
%         end
%         sparse(ML1-full(ML))
    case 'FEM'
        N = pCalc.N; n_mod = pCalc.n_mod; % n_mod=(n-1)*N+1 - number of rows
        matrix = getFEMParams(N);
        b = zeros(1, n-1); c = b; % (n-1) - number of blocks in ML 
        jj=1:(n-1);
        b(jj) = hbar^2/(2*h)*B(pTask, h*(jj-1/2));
        c(jj) = h*(V(pTask, h*(jj-1/2))-1i*2*hbar/tau*rho(pTask, h*(jj-1/2)));
        
        %%
        ML = zeros(n_mod);
        for j=1:(n-1)
            v = ((j-1)*N+1):(j*N+1);
            ML(v, v) = ML(v, v) + c(j)*matrix.C/2 + b(j)*2*matrix.A;
        end
        
%         %% EX22 modification (in case of delta function change)
%         if isfield(pTask, 'V_Q') && isfield(pTask, 'V_x0')
%             x0 = pTask.V_x0;
%             j0 = fix(x0/h);
%             f  = calcFiniteElementValue(2*(x0/h-(j0+1/2)), N+1);
%             v  = ((j0-1)*N+1):(j0*N+1);
%             ML(v, v) = ML(v, v) + pTask.V_Q*(transpose(f(1, :))*f(1, :));
%         end
        
        %%
        MR = conj(ML);
        if pCalc.pLBoundary.bUseDTBC % left boundary
            R = pCalc.pLBoundary.R;
            c_0 = pCalc.pLBoundary.c_0;
            ML(1, 1) = ML(1, 1) - hbar^2*B_inf*c_0*R(1);
            MR(1, 1) = MR(1, 1) - hbar^2*B_inf*c_0*R(2);
        else
            ML(1, 1) = 1; MR(1, 1) = 0;
        end
        if pCalc.pRBoundary.bUseDTBC % right boundary
            R = pCalc.pRBoundary.R;
            c_0 = pCalc.pRBoundary.c_0;
            ML(n_mod, n_mod) = ML(n_mod, n_mod) - hbar^2*B_inf*c_0*R(1);
            MR(n_mod, n_mod) = MR(n_mod, n_mod) - hbar^2*B_inf*c_0*R(2);
        else
            ML(n_mod, n_mod) = 1; MR(n_mod, n_mod) = 0;
        end
end
ML =  sparse(ML);
MR = -sparse(MR);