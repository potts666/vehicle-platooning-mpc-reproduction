%% Sun 2026: corrected ABD simulation and explicit paper-constraint audit
% Requires Optimization Toolbox. Run this script from MATLAB.
% mode='engineering': staggered nonlinear MPC with hard terminal equality,
% bounds and collision constraints. Does NOT claim the paper's string proof.
% mode='paper_audit': literal (14)/(6i), including i=1, for auditing only.
% mode='anchored_string': a reproducible, empirical string-stability
% experiment: follower 1 has a nonzero error envelope; (6i) starts at i=2.
% Table I/II are retained. Unreported scenario/bound choices are marked below.
% No infeasible fmincon output is ever applied; diagnostics are saved first.
% Test overrides: SUN_ABD_TSIM, SUN_ABD_MODE, SUN_ABD_NO_PLOTS.
clear; clc;
assert(exist('fmincon','file')==2,'Optimization Toolbox / fmincon required.');
mode = 'engineering';  % 'engineering', 'paper_audit', or 'anchored_string'.
Tsim = 10;
if ~isempty(getenv('SUN_ABD_MODE')), mode=getenv('SUN_ABD_MODE'); end
if ~isempty(getenv('SUN_ABD_TSIM')), Tsim=str2double(getenv('SUN_ABD_TSIM')); end
assert(any(strcmp(mode,{'engineering','paper_audit','anchored_string'})),'Unknown mode.');
C.literalString=strcmp(mode,'paper_audit');
C.anchoredString=strcmp(mode,'anchored_string');
C.centralInit=~strcmp(mode,'engineering');
C.N=4; C.dt=0.1; C.h=0.05; C.Np=20;
C.d=10; C.g=9.8; C.fR=0.01; C.Cd=0.39;
C.R=[10 5 5]; C.Q=[1 10]; C.alpha=0.6; C.beta=0.2;
tableI=[1183.5 3.0 .35 .57 -245 245;1228.7 3.1 .36 .51 -265 265; ...
        1283.2 3.2 .34 .61 -300 300;1342.8 3.4 .37 .59 -350 350];
C.m=tableI(:,1); C.eta=tableI(:,2); C.r=tableI(:,3);
C.tau=tableI(:,4); C.lo=tableI(:,5); C.hi=tableI(:,6);
% Explicit assumptions: bounds and tPrime not tabulated. The virtual leader
% is at 10 m/s before t=0, jumps to 10.5 m/s at t=0, then decelerates
% uniformly to 9.5 m/s over the next 2 s and remains at that speed.
C.vmin=0; C.vmax=30; C.smin=5; C.smax=15;
C.v0=10; C.vJump=10.5; C.vf=9.5; C.leaderDecelDuration=2;
C.leaderAcceleration=(C.vf-C.vJump)/C.leaderDecelDuration;
C.tPrime=-0.1;
% New, explicit design choice (not a value reported by Sun et al.):
% the first physical follower is anchored to the virtual leader by this
% allowable position-error envelope. The predecessor string starts at i=2.
C.e1max=0.30;
if ~isempty(getenv('SUN_ABD_E1MAX')), C.e1max=str2double(getenv('SUN_ABD_E1MAX')); end
assert(isfinite(C.e1max)&&C.e1max>0,'SUN_ABD_E1MAX must be positive.');
C.tol=2e-6; C.scale=[10;10;350];
C.outdir=fullfile(fileparts(mfilename('fullpath')),'sun_abd_results');
if ~exist(C.outdir,'dir'), mkdir(C.outdir); end
C.opts=optimoptions('fmincon','Algorithm','sqp','Display','off', ...
    'MaxIterations',500,'MaxFunctionEvaluations',80000, ...
    'ConstraintTolerance',1e-8,'OptimalityTolerance',1e-6, ...
    'StepTolerance',1e-10,'ScaleProblem','obj-and-constr');
fprintf('Mode: %s. Literal string constraints: %d; anchored extension: %d\n', ...
    mode,C.literalString,C.anchoredString);
if strcmp(mode,'engineering')
    fprintf('Engineering experiment: (14)/(6i) not enforced; string peaks are measured only.\n');
elseif C.anchoredString
    fprintf(['Anchored-string extension: |e_1| <= %.3f m; (14)/(6i) chain starts at i=2. ' ...
        'The max-error ordering is checked after the run; this is not a claim of the paper theorem.\n'],C.e1max);
end
assert(isfinite(Tsim)&&Tsim>0&&abs(Tsim/C.h-round(Tsim/C.h))<1e-8);
N=C.N; Np=C.Np; dt=C.dt; h=C.h;
Tjump=eq_torque(C.vJump,C); T0=eq_torque(C.v0,C);

%% Initialization
% The engineering warm start is explicitly separate from the literal OCP2.
if C.centralInit
    tc=C.tPrime; hc=C.h; H=round(-tc/hc)+2*(Np+1);
else
    tc=0; hc=C.dt; H=Np+1;
end
xStart=[leader_position(tc,C)-(1:N)*C.d; C.v0*ones(1,N); T0'];
centralObj=@(z) central_cost(z,xStart,tc,hc,H,C);
centralCon=@(z) central_constraints(z,xStart,tc,hc,H,C);
z0=repmat((Tjump./C.hi)',H,1);
[zc,initLog]=solve_checked(centralObj,centralCon,z0(:), ...
    repmat((C.lo./C.hi)',H,1),ones(H,N),C,'initialization');
Uc=reshape(zc,H,N).*C.hi'; Xc=cell(1,N);
for i=1:N, Xc{i}=rollout(xStart(:,i),Uc(:,i),i,hc,C); end
plans=repmat(struct('x',[],'u',[],'t',0),1,N);
if C.centralInit
    % Before either group has performed a local update, make the OCP2 plan
    % available from t=0. This supplies an odd-vehicle reference at t=0.05
    % for the first even-vehicle OCP rather than requesting a future plan.
    zidx=round(-tc/hc)+1;
    for i=1:N
        x(:,i)=Xc{i}(:,zidx); %#ok<SAGROW>
        plans(i).x=Xc{i}(:,zidx:2:zidx+2*Np);
        plans(i).u=Uc(zidx:2:zidx+2*(Np-1),i);
        plans(i).t=0;
    end
    initialUIndex=zidx;
else
    x=xStart;
    for i=1:N
        plans(i).u=Uc(1:Np,i);
        plans(i).x=Xc{i}(:,1:Np+1);
        plans(i).t=0;
    end
    initialUIndex=1;
end

%% Each transmitted PREDICTED plan has its own time stamp.
% A self assumed plan is obtained separately at its next update time.
nEvents=round(Tsim/h); time=(0:nEvents)*h;
xHist=nan(3,N,nEvents+1); xHist(:,:,1)=x;
uHist=nan(N,nEvents); leaderHist=zeros(3,nEvents+1);
leaderHist(1,:)=leader_position(time,C);
leaderHist(2,:)=leader_velocity(time,C);
uHeld=Uc(initialUIndex,:)'; slope=zeros(3,N); updateCount=zeros(1,N);
oddUpdate=false(N,nEvents); evenUpdate=false(N,nEvents);
solverLog=cell(N,nEvents); maxResidual=initLog.violation;
for event=1:nEvents
    t=(event-1)*h;
    % Even first update h; odd first update dt, as in Algorithm 1.
    if event==1, active=[];
    elseif mod(event-1,2)==1, active=2:2:N;
    else, active=1:2:N; end
    snapshot=plans;
    for i=active
        updateCount(i)=updateCount(i)+1;
        own=sample_plan(snapshot(i),t+(0:Np)*dt,i,C);
        ua=sample_inputs(snapshot(i),t+(0:Np-1)*dt,i,C);
        if i==1
            front=desired(0,t+(0:Np)*dt,C);
        else
            front=sample_plan(snapshot(i-1),t+(0:Np)*dt,i-1,C);
        end
        back=[];
        if i<N, back=sample_plan(snapshot(i+1),t+(0:Np)*dt,i+1,C); end
        target=desired(i,t+(0:Np)*dt,C);
        obj=@(z) local_cost(z,x(:,i),own,front,back,target,i,C);
        con=@(z) local_constraints(z,x(:,i),own,front,back,target, ...
            i,t,updateCount(i),C);
        try
            [z,logEntry]=solve_checked(obj,con,ua/C.hi(i), ...
                C.lo(i)/C.hi(i)*ones(Np,1),ones(Np,1),C, ...
                sprintf('vehicle_%d_time_%.2f',i,t));
        catch err
            save(fullfile(C.outdir,'failed_run.mat'),'C','mode','time','xHist', ...
                'uHist','plans','event','i','t','solverLog');
            rethrow(err);
        end
        u=z*C.hi(i); X=rollout(x(:,i),u,i,dt,C);
        % Store the unshifted optimum for the neighbor half-step reference.
        plans(i)=struct('x',X,'u',u,'t',t);
        uHeld(i)=u(1); slope(:,i)=dynamics(x(:,i),uHeld(i),i,C);
        solverLog{i,event}=logEntry; maxResidual=max(maxResidual,logEntry.violation);
        if mod(i,2), oddUpdate(i,event)=true; else, evenUpdate(i,event)=true; end
    end
    for i=1:N
        first=h*(1+mod(i,2));
        if t<first-1e-9
            if C.centralInit
                idx=round((t-tc)/hc)+1;
                uHeld(i)=Uc(idx,i);
                slope(:,i)=dynamics(x(:,i),uHeld(i),i,C);
            else
                if event==1, slope(:,i)=dynamics(x(:,i),uHeld(i),i,C); end
                x(:,i)=x(:,i)+h*slope(:,i);
            end
        else
            % Partial Euler segment uses frozen slope from the vehicle's
            % update instant: two half segments equal the paper dt Euler step.
            x(:,i)=x(:,i)+h*slope(:,i);
        end
    end
    xHist(:,:,event+1)=x; uHist(:,event)=uHeld;
    actualGap=[leader_position(t+h,C),x(1,1:end-1)]-x(1,:);
    physicalViolation=max([0, C.smin-actualGap, actualGap-C.smax, ...
        C.vmin-x(2,:),x(2,:)-C.vmax,(C.lo'-x(3,:))/350,(x(3,:)-C.hi')/350]);
    if physicalViolation>1e-4
        save(fullfile(C.outdir,'physical_violation.mat'));
        error('Actual bounds violated at %.2f s (%.3g); run stopped.',t+h,physicalViolation);
    end
    if mod(event,20)==0, fprintf('t=%.1f s, max scaled residual %.2g\n',t+h,maxResidual); end
end

%% Measurements are not mathematical proofs.
positions=squeeze(xHist(1,:,:));
trackingError=positions-leaderHist(1,:)+(1:N)'*C.d;
gap=[leaderHist(1,:);positions(1:end-1,:)]-positions;
peak=max(abs(trackingError),[],2);
stringPass=all(peak(2:end)<=peak(1:end-1)+1e-5);
fprintf('Peak error [m]: %s\n',mat2str(peak',5));
fprintf('Final error [m]: %s\n',mat2str(trackingError(:,end)',5));
fprintf('Gap min/max [m]: %.6f / %.6f\n',min(gap,[],'all'),max(gap,[],'all'));
fprintf('Measured follower peak ordering: %d (not a proof).\n',stringPass);
fprintf('Max accepted scaled constraint residual: %.3g\n',maxResidual);
save(fullfile(C.outdir,[mode '_results.mat']));
if isempty(getenv('SUN_ABD_NO_PLOTS'))
    fig1=figure('Color','w','Name',['Sun ABD: ' mode]);
    subplot(3,1,1); plot(time,trackingError','LineWidth',1.2); grid on
    ylabel('tracking error [m]'); title(['Mode: ' mode],'Interpreter','none');
    legend(compose('vehicle %d',1:N),'Location','eastoutside');
    subplot(3,1,2); plot(time,[leaderHist(2,:);squeeze(xHist(2,:,:))]','LineWidth',1.2); grid on
    ylabel('speed [m/s]'); legend(["leader",compose('vehicle %d',1:N)]);
    subplot(3,1,3); stairs(time(1:end-1),uHist','LineWidth',1.1); grid on
    ylabel('command torque [N m]'); xlabel('time [s]');
    fig2=figure('Color','w','Name','Gap and schedule');
    subplot(2,1,1); plot(time,gap'); hold on; yline(C.smin,'r--'); yline(C.smax,'r--');
    yline(C.d,'k:'); grid on; ylabel('gap [m]');
    subplot(2,1,2); stairs(time(1:end-1),[oddUpdate(1,:);evenUpdate(2,:)]');
    grid on; xlim([0 min(1,Tsim)]); xlabel('time [s]'); legend('odd','even');
    savefig(fig1,fullfile(C.outdir,[mode '_tracking.fig']));
    savefig(fig2,fullfile(C.outdir,[mode '_gaps.fig']));
    exportgraphics(fig1,fullfile(C.outdir,[mode '_tracking.png']),'Resolution',150);
    exportgraphics(fig2,fullfile(C.outdir,[mode '_gaps.png']),'Resolution',150);
end

function [z,info]=solve_checked(obj,con,z0,lb,ub,C,label)
lb=lb(:); ub=ub(:); z0=z0(:);
[c0,e0]=con(z0);
if max([0;c0(:);abs(e0(:));lb-z0;z0-ub])>C.tol
    % Restore nonlinear feasibility before minimizing the nonsmooth norm
    % objective. This does not alter any constraint or add a slack variable.
    restoreOpts=optimoptions('lsqnonlin','Display','off','MaxIterations',300, ...
        'MaxFunctionEvaluations',60000,'FunctionTolerance',1e-14, ...
        'StepTolerance',1e-12,'OptimalityTolerance',1e-12);
    z0=lsqnonlin(@(v) feasibility_residual(v,con),z0,lb,ub,restoreOpts);
end
[c0,e0]=con(z0); seedViolation=max([0;c0(:);abs(e0(:));lb-z0;z0-ub]);
if seedViolation>C.tol
    info=struct('violation',seedViolation); c=c0; e=e0; z=z0;
    save(fullfile(C.outdir,'solver_failure.mat'),'z','z0','c','e','info','label','C');
    error('SUN:InitializationInfeasible','%s: feasibility restoration failed (%.3g). No input applied.',label,seedViolation);
end
[z,f,flag,out]=fmincon(obj,z0,[],[],[],[],lb,ub,con,C.opts);
[c,e]=con(z); violation=max([0;c(:);abs(e(:));lb-z;z-ub]);
info=struct('exitflag',flag,'violation',violation,'objective',f,'output',out,'fallback',false);
if flag<=0 || violation>C.tol || ~all(isfinite(z)) || ~isfinite(f)
    % Only the independently checked feasible seed can be a fallback.
    z=z0; f=obj(z); [c,e]=con(z);
    violation=max([0;c(:);abs(e(:));lb-z;z-ub]);
    info.fallback=true; info.violation=violation; info.objective=f;
    warning('SUN:FeasibleFallback','%s: optimizer exit %d; applying verified feasible candidate, not an optimum.',label,flag);
end
if ~all(isfinite(z)) || ~isfinite(f) || ~isfinite(violation) || violation>C.tol
    save(fullfile(C.outdir,'solver_failure.mat'),'z','z0','c','e','info','label','C');
    error('SUN:SolverFailure','%s: exitflag=%d, scaled violation=%.3g. No input applied. %s', ...
        label,flag,violation,out.message);
end
end

function r=feasibility_residual(z,con)
[c,e]=con(z); r=[max(c(:),0);e(:)];
end

function J=local_cost(z,x,xa,xf,xb,xd,i,C)
X=rollout(x,z*C.hi(i),i,C.dt,C); J=0;
assert(isequal(size(X),size(xd)),'Desired trajectory dimension mismatch.');
% Delta x includes equilibrium torque differences for heterogeneous vehicles.
tf=desired(i-1,0,C)-desired(i,0,C);
for j=1:C.Np
    J=J+C.R(1)*norm(X(:,j)-xa(:,j))+C.R(2)*norm(X(:,j)-xf(:,j)+tf);
    if ~isempty(xb)
        tb=desired(i,0,C)-desired(i+1,0,C);
        J=J+C.R(3)*norm(xb(:,j)-X(:,j)+tb);
    end
end
% Hard terminal equality is imposed below. No matrix-valued terminal error.
end

function [c,e]=local_constraints(z,x,xa,xf,xb,xd,i,t,count,C)
X=rollout(x,z*C.hi(i),i,C.dt,C);
c=state_bounds(X,i,C);
gap=xf(1,2:end)-X(1,2:end);
c=[c;(C.smin-gap)'/10;(gap-C.smax)'/10];
if ~isempty(xb)
    gap=X(1,2:end)-xb(1,2:end);
    c=[c;(C.smin-gap)'/10;(gap-C.smax)'/10];
end
if C.literalString
    epredFront=xf(1,1:C.Np)-leader_position(t+(0:C.Np-1)*C.dt,C)+(i-1)*C.d;
    c=[c;(abs(X(1,1:C.Np)-xa(1,1:C.Np))-C.beta^count*abs(epredFront))'/10];
elseif C.anchoredString
    % The virtual leader is a reference, not a predecessor carrying a
    % nonzero tracking error. Anchor follower 1 explicitly; for i>=2 retain
    % the paper's predicted-vs-assumed mismatch condition (6i).
    if i==1
        e1=X(1,1:C.Np)-leader_position(t+(0:C.Np-1)*C.dt,C)+C.d;
        c=[c;(abs(e1)-C.e1max)'/10; ...
            (abs(X(1,1:C.Np)-xa(1,1:C.Np))-C.beta^count*C.e1max)'/10];
    else
        epredFront=xf(1,1:C.Np)-leader_position(t+(0:C.Np-1)*C.dt,C)+(i-1)*C.d;
        c=[c;(abs(X(1,1:C.Np)-xa(1,1:C.Np))-C.beta^count*abs(epredFront))'/10];
    end
end
e=(X(:,end)-xd(:,end))./C.scale;
end

function J=central_cost(z,x,t,h,H,C)
U=reshape(z,H,C.N).*C.hi'; J=0;
te=eq_torque(leader_velocity(t+(0:H-1)*h,C),C);
for i=1:C.N
    X=rollout(x(:,i),U(:,i),i,h,C); xd=desired(i,t+(0:H)*h,C);
    J=J+C.Q(1)*sum(abs(U(:,i)-te(i,:).'))+C.Q(2)*sum(vecnorm(X(:,1:H)-xd(:,1:H)));
end
end

function [c,e]=central_constraints(z,x,t,h,H,C)
U=reshape(z,H,C.N).*C.hi'; c=[]; e=[]; prev=desired(0,t+(0:H)*h,C);
prevError=zeros(1,H+1);
for i=1:C.N
    X=rollout(x(:,i),U(:,i),i,h,C); xd=desired(i,t+(0:H)*h,C);
    gap=prev(1,2:end)-X(1,2:end);
    c=[c;state_bounds(X,i,C);(C.smin-gap)'/10;(gap-C.smax)'/10]; %#ok<AGROW>
    terminal=(X(:,end-1:end)-xd(:,end-1:end))./C.scale;
    % At constant target speed, terminal p/v at the last point follow
    % exactly from the preceding equilibrium state and Euler dynamics.
    % Keep the 4 independent equalities, avoiding a rank-deficient SQP QP.
    e=[e;terminal(:,1);terminal(3,2)]; %#ok<AGROW>
    err=X(1,:)-xd(1,:);
    if C.literalString
        % Literal reading of (14), including i=1 against virtual CAV 0.
        c=[c;(abs(err)-C.alpha*abs(prevError))'/10]; %#ok<AGROW>
    elseif C.anchoredString
        if i==1
            % New boundary condition: a finite error allowance for the first
            % physical follower. It replaces the ill-posed e_1 <= alpha e_0.
            c=[c;(abs(err)-C.e1max)'/10]; %#ok<AGROW>
        end
        % Do not falsely replace the max-over-time definition (4c) by a
        % pointwise hard condition for i>=2. That stronger condition is kept
        % only in paper_audit. Actual peak ordering is reported post-run.
    end
    prev=X; prevError=err;
end
end

function c=state_bounds(X,i,C)
Y=X(:,2:end);
c=[(C.vmin-Y(2,:))'/10;(Y(2,:)-C.vmax)'/10; ...
    (C.lo(i)-Y(3,:))'/350;(Y(3,:)-C.hi(i))'/350];
end

function Y=sample_plan(plan,times,i,C)
Y=zeros(3,numel(times));
for j=1:numel(times)
    a=(times(j)-plan.t)/C.dt;
    assert(a>=-1e-7,'Reference requested before its timestamp.');
    k=min(floor(a+1e-8),size(plan.x,2)-1); rem=times(j)-(plan.t+k*C.dt);
    xx=plan.x(:,k+1);
    if k<numel(plan.u)
        u=plan.u(k+1);
    else
        te=eq_torque(leader_velocity(times(j),C),C); u=te(i);
    end
    Y(:,j)=xx+rem*dynamics(xx,u,i,C);
end
end

function u=sample_inputs(plan,times,i,C)
te=eq_torque(leader_velocity(times,C),C); u=te(i,:).';
for j=1:numel(times)
    k=floor((times(j)-plan.t)/C.dt+1e-8)+1;
    if k>=1&&k<=numel(plan.u), u(j)=plan.u(k); end
end
end

function X=rollout(x,u,i,h,C)
X=zeros(3,numel(u)+1); X(:,1)=x;
for j=1:numel(u), X(:,j+1)=X(:,j)+h*dynamics(X(:,j),u(j),i,C); end
end
function f=dynamics(x,u,i,C)
f=[x(2);(C.eta(i)/C.r(i)*x(3)-C.Cd*x(2)^2)/C.m(i)-C.g*C.fR;(u-x(3))/C.tau(i)];
end
function T=eq_torque(v,C)
T=C.r./C.eta.*(C.Cd*v^2+C.m*C.g*C.fR);
end
function xd=desired(i,t,C)
v=leader_velocity(t,C);
if i==0
    torque=zeros(size(t));
else
    te=eq_torque(v,C); torque=te(i,:);
end
xd=[leader_position(t,C)-i*C.d;v;torque];
end
function v=leader_velocity(t,C)
v=C.vf*ones(size(t));
past=t<0; decelerating=t>=0 & t<C.leaderDecelDuration;
v(past)=C.v0;
v(decelerating)=C.vJump+C.leaderAcceleration*t(decelerating);
end
function p=leader_position(t,C)
p=zeros(size(t));
past=t<0; decelerating=t>=0 & t<=C.leaderDecelDuration;
p(past)=C.v0*t(past);
p(decelerating)=C.vJump*t(decelerating) ...
    + 0.5*C.leaderAcceleration*t(decelerating).^2;
pAtEnd=C.vJump*C.leaderDecelDuration ...
    + 0.5*C.leaderAcceleration*C.leaderDecelDuration^2;
steady=t>C.leaderDecelDuration;
p(steady)=pAtEnd+C.vf*(t(steady)-C.leaderDecelDuration);
end
