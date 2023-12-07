clear all
close all
clc

%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%
omega=10; %FREQUENCY HO
ga=0;omega; gb=0;omega*1.4; %COUPLINGS
number_of_oscillations=1.5; trotter_steps=300; %SIMULATION PARAMS
omega_a=4.403989664625457e-05; omega_b=48*omega_a;%SPIN FREQUENCIES

%Simulation time 
tf=number_of_oscillations*2*pi/omega; %final time
dt=tf/trotter_steps; %required time step
n=round(tf/dt); %number of times we will apply the loop to get to tf 
tt=linspace(0,tf,n); %timevector
dt=tt(2)-tt(1); %actual timestep
Plotting=1; %Plottiong the dynamics 
wigner=0;  %Plottiong the dynamics with wigner 

%%%%%%%%%%%%%%%% Definition of the grid  %%%%%%%%%%%%%%
nx=2^7; Lx=(nx-1)*(2*pi/nx)^.5; 
Lx0=-Lx/2; 
dx=Lx/(nx-1); x=Lx0+dx*(0:nx-1); dpx=(2*pi)/(nx*dx); 
px=dpx*(0:nx/2);px_aux=[-px(end:-1:2),px(1:end-1),];
px=[px(1:end-1),-px(end:-1:2)];


%%%%%%%%%%%%%%%%%%% Quantities tracked %%%%%%%%%%%%%%%%%%
P_spin_a=[]; P_spin_b=[]; P_spin_ho=[]; %PURITIES
LNs=[]; %LOGARITHMIC NEGATIVITY 


%%%%%%%%%%%%%%%%%%%%%% INITIAL STATE %%%%%%%%%%%%%%%%%%%%%%

%Coherent sate for x
psi0_qho=exp(-(x-Lx/7).^2/(2))/norm(exp(-(x-Lx/7).^2/(2))); 

%Spin states: superposition
psi0_spin_1=1/2^.5*[1,1i];
psi0_spin_2=1/2^.5*[1,1i];

%Total state
psi0=kron(psi0_spin_2,kron(psi0_spin_1,psi0_qho));


%%%%%%%%%%%%%%%%% Hamiltonian operators %%%%%%%%%%%%%%%%
spin_a=[ones(1,2*nx),-ones(1,2*nx)];
spin_b=[ones(1,nx),-ones(1,nx),ones(1,nx),-ones(1,nx)];
spin_a_qho=-(2)^.5*ga*[x.*ones(1,nx),x.*ones(1,nx),-x.*ones(1,nx),-x.*ones(1,nx)];
spin_b_qho=-(2)^.5*gb*[x.*ones(1,nx),-x.*ones(1,nx),x.*ones(1,nx),-x.*ones(1,nx)];
kin=[omega/2*px.^2,omega/2*px.^2,omega/2*px.^2,omega/2*px.^2]; %kinetic vector
V=[1/2*omega*x.^2,1/2*omega*x.^2,1/2*omega*x.^2,1/2*omega*x.^2]; %potential vector



%%%%%%%%%%%%%%%% Initial record of quantities %%%%%%%%%%%%%%%%
psi=psi0;
psi_aux=psi.';
rho= psi_aux*psi_aux';

%Purities
XPT_spin_a = PartialTrace(rho,[2,3],[2,2,nx],-1);
XPT_spin_b = PartialTrace(rho,[1,3],[2,2,nx],-1);
XPT_ho = PartialTrace(rho,[1,2],[2,2,nx],-1);
P_spin_a(1)=trace(XPT_spin_a* XPT_spin_a);
P_spin_b(1)=trace(XPT_spin_b* XPT_spin_b);
P_spin_ho(1)=trace(XPT_ho* XPT_ho);

%Logarithmic negativity between spins
rho_12 = PartialTrace(rho,2,[4,nx],-1);
rho_T2 = PartialTranspose(rho_12,2,[2,2]);
LN=log2(sum(eig(rho_T2*rho_T2').^.5));
LN=max([real(LN),0]);
LNs(1)=LN;


%%%%%%%%%%%%%%%% Time Evolution %%%%%%%%%%%%%%%%
for i=1:n-1

    disp(i)
  

    psi=psi.*exp(-1i*V*dt).*exp(-1i*omega_a*spin_a*dt).*exp(-1i*omega_b*spin_b*dt).*exp(-1i*spin_a_qho*dt).*exp(-1i*spin_b_qho*dt);
    psip=[fft(psi(1:nx)),fft(psi(nx+1:2*nx)),fft(psi(2*nx+1:3*nx)),fft(psi(3*nx+1:end))];
    psip=psip.*exp(-1i*dt*kin); 
    psi=[ifft(psip(1:nx)),ifft(psip(nx+1:2*nx)),ifft(psip(2*nx+1:3*nx)),ifft(psip(3*nx+1:end))]; %inversetransform 
    

    if Plotting==1

     plot(x,abs(psi(1:nx)).^2+abs(psi(nx+1:2*nx)).^2+abs(psi(2*nx+1:3*nx)).^2+abs(psi(3*nx+1:end)).^2)
     %ylim([0 0.15])
     drawnow

    end 
   

    if wigner==1
        psi_aux=psi.';
        rho=psi_aux* psi_aux';
        for j=1:4
            for l=1:4
                figure(3)
                subplot(4,4,l+4*(j-1))
                sub_rho= rho(1+(j-1)*nx:j*nx, 1+(l-1)*nx:l*nx);
                result=RhoWigner( sub_rho,x,px_aux);
                imagesc(real(result))
                %colorbar 
                figure(4)
                subplot(4,4,l+4*(j-1))
                imagesc(imag(result))
                %colorbar 
            end 
        end 
        figure(3)
        cm = 'Purples4';
        colormap(othercolor(cm))
        colormap(flipud(othercolor(cm)))
        sgtitle('Real part','Interpreter','latex','FontSize',30) 
        figure(4)
        cm = 'Purples4';
        colormap(othercolor(cm))
        colormap(flipud(othercolor(cm)))
        sgtitle('Imaginary part','Interpreter','latex','FontSize',30) 
            
   end 





    %% record of quantities %%
    psi_aux=psi.';
    rho= psi_aux*psi_aux';

    XPT_spin_a = PartialTrace(rho,[2,3],[2,2,nx],-1);
    XPT_spin_b = PartialTrace(rho,[1,3],[2,2,nx],-1);
    XPT_ho = PartialTrace(rho,[1,2],[2,2,nx],-1);
    P_spin_a(i+1)=trace(XPT_spin_a* XPT_spin_a);
    P_spin_b(i+1)=trace(XPT_spin_b* XPT_spin_b);
    P_spin_ho(i+1)=trace(XPT_ho* XPT_ho);

    rho_12 = PartialTrace(rho,2,[4,nx],-1);
    rho_T2 = PartialTranspose(rho_12,2,[2,2]);
    LN=log2(sum(eig(rho_T2*rho_T2').^.5));
    LN=max([real(LN),0]);
    LNs(i+1)=LN;

   
end

%%%%%%%%%%%%%%% plotting retuls %%%%%%%%%%%%%
colors=[[0.97,0.91,0.99];[0.85,0.76,0.89];[0.67,0.53,0.73];[0.51,0.33,0.59];[0.6836,    0.6641,    0.6914]; [0.44,0.23,0.53];[0.30,0.09,0.39]];
figure(2)
hold on

plot(linspace(0,1.5,trotter_steps),P_spin_a, '--','color',colors(7,:),'LineWidth',1.5)
plot(linspace(0,1.5,trotter_steps),P_spin_b, '-','color',colors(5,:),'LineWidth',1.5)
plot(linspace(0,1.5,trotter_steps),P_spin_ho, '-.','color',colors(3,:),'LineWidth',1.5)
legend({'Purity_spin_a','Purity_spin_b','Purity_ho'},'Interpreter','latex','FontSize',10)
xlabel('$t/\omega$', 'FontSize',20,'Interpreter','latex')
ylabel('Purity', 'FontSize',20,'Interpreter','latex','Color','k')
legend boxoff

figure(3)
plot(linspace(0,1.5,trotter_steps),LNs,'-','color',colors(3,:),'LineWidth',1.5)
xlabel('$t/\omega$', 'FontSize',20,'Interpreter','latex')
ylabel('Logarithmic Negativity', 'FontSize',20,'Interpreter','latex','Color','k')
legend({'Logarithmic negativity'},'Interpreter','latex','FontSize',20)
hold on
ax = gca;
%ylim([-0.1,0.7])
legend boxoff


