% Comments after main
% file name is bliss20yya-z.m based on year and version a-z.
% version 2019b started 2019/10/31
% version 2020a started 2020/01/17
% requires directory graphname/ with files 
% version 2020b started 2020/07/02
% Morse Index list added
% version 2023a started 2023/06/02
% branch colors changed so trivial solution (j = 1) and invariant subspace
% = R^N (j = niss) are black, and the others cycle through colors
% version 2023b started 2023/08/07
% The Laplacian (lap) is changed to the Weighted Adjacency Matrix (wam)

% aut.txt List of automorphisms (permutations of 1 to N) including
% identity.

% wam.txt the weighted adjacency matrix in sparse format, 
% i j M_{i,j}.  If M = 0, a file
% with one line 1 1 0 is allowed.

% inv.txt with the invariant subspaces for odd f. The j'th row is the color
% vector for W_j.  The first row must be 0 0 ... 0, and it is best to list
% the subspaces with non-decreasing dimension.  The strict requirement is
% if W_j \subsetneq W_k then
% colorlist for W_j should preceed colorlist for W_k in inv.txt.)

% bif.txt with subsets:  line m  * * d * * means W_m \subseteq W_d

function BLIS_2023b
clear;
global jobs bcnt speed smin smax verbosity blist
global maxLevel maxLevelKickedIn ymax
global wam f fp dfds MImax f_is_odd y graphName auts print_plist

% frequently changed parameters
graphName = ['diamond'];

SaveThisRun = false; % saves .txt of Command Window and .pdf of figure.
                     % Answer "y" in command window to really save files.

speed = .1;     % distance between computed points in bifurcation diagram

maxLevel = 1; % Follow daughters to this level. (First branch is level 0.)

plot_bifpoints = true;
plot_foldpoints = false;
plot_branch_as_line = true; % false useful for setting speed

smin = -4; smax = 4.5;        % bifurcation parameter window
ymax = 2;  %default 1000. Needed with asymptotically linear nonlinearity

maxNumBranches  = 50;  %  Default is 50

verbosity       = 1;       % output control: 0 is minimal, 1 gives more
                           % info about daughters, 2 gives lots of info
print_plist     = false;   % print s y MI newton_its at end of branch

startWithTrivialBranch = true;
if startWithTrivialBranch
    sInit = smin;  % default
    %sInit = 70;  % un-comment this line to set sInit != smin
    usInit = [sInit]; vsInit = [1]; jInit = 1; % trivial branch
    % Needs sInit >= smin with vsInit = [1] (going right).
    % or sInit <= smax with vsInit = [-1] (going left).
    
else
    %usInit = [0.9003; -1.63; 2]; vsInit = [-1; -1; 3]; jInit = 6; % diamond
    usInit = [1.732; -3]; vsInit = [-1; 3]; jInit = 2; % diamond
    
    % must be used if f(0,s) != 0.  Start with any approx or exact solution
    % in the invariant subspace W_jInit with dimension d 
    % usInit is the (approximate) initial solution [u1; ... ud; s]
    % vsInit is the (approximate) tangent vector [v1; ... vd; s] pointing
    % into the plotting region
end

% The following parameters are more spread out

% You can change the label for the schematic function, using LaTex.  
% If you do, search for 'switch y' and change those y strings.
% The value of those y strings needs to agree with the y defined here.

% CHOOSE SCHEMATIC FUNCTION WITH THE NEXT LINE
switch 4 % change integer to get desired schematic function y vs. s
    % don't change the case statements unless you know what you're doing
    % Many schematice functions put two branches on top of each other.
    % Choice "switch 6" makes sure branches don't lie on top of each other. 

    % Don't change the numbers below this line!
    case 1 % This defines y in case the line above is "switch 1"
        y = 'u_1';

    case 2
        y = 'u_N';

    case 3
        y = '\frac 1 N \|u\|_2^2';

    case 4
        y = '\frac 1 N \|u\|_1';

    case 5
        y = '\frac 1 {\sqrt N} \|u\|_2';

    case 6   % this one separates branches.
        y = '\displaystyle{\sum_{i=1}^N \sin(i) x_i}';  

    case 7
        y = '\|u\|_\infty'; 
        
    case 8
        y = 'N u_N';


end

% the derivatives fp = D_1 f and dfds = D_2 f need to be input by hand.
% Usually f(u, s) = s*u + nonlinearity(u)

% CHOOSE NONLINEARITY WITH THE NEXT LINE
switch 1 % Change integer and possibly c5 or c2 to get desired nonlinearity
    
    case 1   % f is odd, this is for "switch 1" in the previous line
        c5 = 0; % coefficient of degree 5 term.
        % Constant branch has fold at s = -1/4 * 1/c5 = -1/(4 c5)
        f    = @(u,s) s*u + u.^3 - c5 * u.^5;
        fp   = @(u,s) s + 3*u.^2 - 5 * c5 * u.^4;
        dfds = @(u,s) u;
        
    case 2  % f is not odd if c2 ~= 0
        c2 = 0; % coef of degree 2 term. 
        % Constant branch has fold at s = -c2^2/4
        f    = @(u,s) s*u - u.^3 + c2*u.^2;
        fp   = @(u,s) s - 3*u.^2 + 2*c2*u;
        dfds = @(u,s) u;
        
    case 3   % f(0, s) ~= 0
        f    = @(u,s) s*u + u.^3 + 1.;
        fp   = @(u,s) s + 3*u.^2;
        dfds = @(u,s) u;
        
    case 4   % df/ds(0,s) ~= u
        f    = @(u,s) (s*s-1)*u - u.^3;% - 1/2*u.^5;
        fp   = @(u,s) (s*s-1) - 3*u.^2;%  - 5/2* u.^4;
        dfds = @(u,s) 2*s*u;
        
    case 5    % f us asymptotically linear
        cs = -1;  % This gives vertical asymptote at s = lambda + cs
        f    = @(u,s) s*u + cs*(tanh(u) - u);
        fp   = @(u,s) s + cs*(sech(u).^2 -1);
        dfds = @(u,s) u;

    case 6    % gives non-BLIS bifurcation
        f    = @(u,s) s*s - u.^2;
        fp   = @(u,s) -2*u;
        dfds = @(u,s) 2*s;
end


%%%% end of parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f_is_odd = (f(1.125,0)+f(-1.125,0)==0)&&(f(2.5, .25)+f(-2.5,.25) == 0);
% This can be over-written by true or false by including a line like
% f_is_odd = true;
% The input of 1.125 is chosen since it
% terminates in base 2, and because it is very unlikely to accidentally
% give a true.  (Choosing an integer might cause a problem.)

%work_file       = 'temp';       % for branchmanager, not currently used

maxLevelKickedIn = false;

if f_is_odd
    fOddStr = ', f odd';
else
    fOddStr = ', f not odd';
end
% 2021-04-30
% I might want to have 3 levels, 'f odd', 'f(0) not 0', 'f(0) = 0, f not
% odd'  These could be integers.  We might want to set this with the case
% statement (switch) that chooses f.
% We might be able to use these 3 levels to change the bif.txt information,
% and keep the W_j nunbers the same  (currently the W_j determines the
% color of the branch.  Handling f(0) not 0 might be difficult.

%%%% start output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if SaveThisRun
    prompt = 'Do you really want to save this run? y/n [n]: ';
    str = input(prompt,'s');
    if isempty(str)
        str = 'n';
    end
    if str == 'y'
        fprintf('OK, saving run.\n');
    else
        SaveThisRun = false;
        fprintf('initialize SaveThisRun with the other parameters.\n');
    end
end
tic; % start clock for elapsed time printed with toc 
    
DateStr = datestr(now, 'yyyy-mm-dd_HH-MM');
%This saves a "diary", a copy of the command window to a dated file in the
%directory with the graph name
if SaveThisRun
    fprintf('------------------------------------------------------ ');
    fprintf('Run will be saved to a file since SaveThisRun is true\n');
    OutFileName = [graphName '/' DateStr];
    diary ([OutFileName '.txt']);
else
    fprintf('------------------------------------------------------ ');
    fprintf('Run will not be saved, since SaveThisRun is false\n');
end
figTitle = [graphName fOddStr];

fprintf('-------------------------------------------------------------\n');
fprintf(['Solve f(u_i,s) - (M u)_i = 0, for i = 1, ..., N, \n']);
fprintf(['where M is a Weighted Adjacency Matrix \n']);
fprintf(['The unknown u is in R^N , and the function f: R^2 -> R is \n'])
fprintf(['threaded over u to get a function f: R^N x R -> R^N. \n'])
fprintf(['So the program really solves f(u,s) - Mu = 0.\n\n'])
fprintf(['This run was made at date_time ' DateStr '\n']);
fprintf(['The 4 input .txt files are in the directory ' graphName '\n'])
fprintf(['f(u,s) = '])
fstructure = functions(f);
display(fstructure.function);
if exist('c5', 'var')
    fprintf('c5 = %g\n', c5);
end
if exist('c2', 'var')
    fprintf('c2 = %g\n', c2);
end
if exist('cs', 'var')
    fprintf('cs = %g\n', cs);
end


fprintf('smin = %.5g,  smax = %.5g, speed = %.5g\n', ...
    smin, smax, speed);
fprintf('maxNumBranches = %d,  maxLevel = %d, verbosity = %d.\n', ...
    maxNumBranches, maxLevel, verbosity);

if length(usInit) == 1
    fprintf('Start on trivial branch at s = %.5g with vs = %.5g\n', ...
        usInit(1), vsInit(1));
else
    fprintf('usInit = [ ');
    fprintf('%.5g ', usInit(1:end));
    fprintf('], vsInit = [ ');
    fprintf('%.5g ', vsInit(1:end));
    fprintf('], jInit = %d', jInit);
end

%auts = readmatrix([graphName '/aut.txt'])
auts = importdata([graphName '/aut.txt']);
[~,N] = size(auts);

fprintf('\nN = %d\n', N);

sps = importdata([graphName '/wam.txt']);
wam = sparse(sps(:,1),sps(:,2),sps(:,3), N, N);
 
MImax = min(20, N);  % This is hardwirwed but we might want to change it.

% The structure iss contains information about Invariant SubSpaces
make_iss;

lam = evals([0], 1); % The eigenvalues of the trivial solution at s = 0.
if(abs(lam(1)) < .1^15)
    lam(1) = 0;  % This makes output prettier if first eval is near 0
end
fprintf('The first %i eigenvalues of the Weighted Adjacency Matrix are\n' ...
    , MImax);
if(norm(imag(lam(:))) < .1^10)
    fprintf('%.4g   ', lam(:));
else
    fprintf('%.4g%+.4gi   ', [real(lam(:)), imag(lam(:))].');
    fprintf('\nThe real parts of the first %i eigenvalues are \n', MImax);
    fprintf('%.4g   ', lam(:));
end
fprintf('\n');
%%%%%%%%%%%%%%

% I don't understand why both these lines are needed
jobs = struct('us',usInit,'vs',vsInit,'lev',0,'hist',[], 'issj', jInit);
jobs(1) = struct('us',usInit,'vs',vsInit,'lev',0,'hist',[], 'issj', jInit);

blist = struct('startpt',[],'bifpts',[],'lev',[],'pts',[],'mi',[],...
    'lastnewt',[],'lastsec',[],'hist',[],'lastpt',[], 'issj', [],'foldpts',[]);
bcnt = 0;


while ~isempty(jobs) && bcnt < maxNumBranches
    bcnt = bcnt + 1;
    us = jobs(1).us; vs = jobs(1).vs; lev = jobs(1).lev; ...
        hist = jobs(1).hist; issj = jobs(1).issj; jobs(1) = [];
    fprintf('\n-------- branch %d -------\n',...
        bcnt );
    sInit = full(us(end));
    if abs(sInit) < .1^10
        sInit = 0;
    end
    fprintf('Starting branch %d at s = %.5g, ',...
        bcnt, sInit );
        %bcnt, full(us(end)) );

    fprintf('y = %.5g, in invariant subspace W_%d', ...
        schematic(full(us), issj), issj);
    fprintf([' (',colorString(issj),').\n']);
    newHist = [hist, [bcnt; issj]];
    fprintf('Array shows history of branch numbers, ');
    fprintf('with invariant subspace below:\n');
    disp(newHist);
    blist(bcnt).startpt = us;
    blist(bcnt).lev     = lev;
    blist(bcnt).hist    = newHist;    
    
    branch(us, vs, issj);
    %     branch_old(w,v,ju,gu);
     
    if length(blist(bcnt).pts) >= 2
        % plot branch
        figure(1);
        if plot_branch_as_line
            plot(blist(bcnt).pts(:,1),blist(bcnt).pts(:,2), ...
                'Color', colorRGB(issj), 'LineWidth', 1.5);
        else % plot computed points as dots
            plot(blist(bcnt).pts(:,1),blist(bcnt).pts(:,2), ...
                '.', 'Color', colorRGB(issj));
        end
        xlim([smin smax]); %ylim([0 13]); % JWS took out ylim
        title(figTitle, 'FontSize', 20);
        xlabel('$s$', 'Interpreter', 'latex', 'FontSize', 20); % was 14
        ylabel(['$' y '$'], 'Interpreter', 'latex', 'FontSize', 20);
        hold on;
    else
        % delete branch from structure
        fprintf('***** deleting empty branch %d *****\n',bcnt);
        blist(bcnt) = [];
        bcnt = bcnt - 1;
    end
    
end
if bcnt >= maxNumBranches
    fprintf('\nmaxNumBranches = %d has been reached.\n', maxNumBranches);
end

if maxLevelKickedIn
    fprintf('\nIncrease maxLevel from %d to follow more daughters.\n', ...
        maxLevel);
end
if plot_bifpoints
    for i = 1:length(blist)
        bp = blist(i).bifpts;
        if length(bp) > 1
            plot(bp(:,1),bp(:,2),'o', 'Color', [.6, .6, .6]); % o is circle
        end
    end
end

if plot_foldpoints
    for i = 1:length(blist)
        fpts = blist(i).foldpts;
        if length(fpts) > 1
            plot(fpts(:,1),fpts(:,2),'+', 'Color', [.4,.4,.4]); % + is cross
        end
    end
end

hold off;

if SaveThisRun
    %saveas(figure(1), [OutFileName '.pdf']);
    orient(figure(1),'landscape');
    print(figure(1), [OutFileName '.pdf'], '-dpdf', '-fillpage');
    %print(figure(1), [OutFileName '.pdf'], '-dpdf');
    diary off;
end
%figure(3); plotuT(u,T); % my so-so contour plot
%figure(3); plotu(iss(issj).B*us(1:end-1), T, verts);

%save(work_file);
toc % elapsed time for the whole program.
end %%%%%%%%%%%%%%%%%%%%%%%%% of main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% versions blis2020b.m started 2020/08/02. The new version is due to the
    % introduction of the MI_list.
% version blis2020a.m 
% 2020/07/17. A few changes over the last few days.
    % in make_iss, the Weighted Adjacency Matrix iss(j).wam is no longer rounded.  
    % This caused errors when I used a scaled Weighted Adjacency Matrix, 
    % which has floating point dx^2 prefactor.   
    % iss(j).B, the {0, 1, -1} basis matrix has a bug when size oa A_i and
    % A_{-i} are not the same.  This has been fixed
%
% version "b", blis2019b.m
% The bif.txt is no longer signed to
% indicate pitchforks, but the computation isPitchfork is used.  The last a
% version had both methods and compared them.  I have saved "bifSigned.txt"
% for use with that old matlab program.  Note how I saved bifOdd.txt and
% bifNotOdd.txt in tadpole/.  I would move the correct one into bif.txt and
% then run it, but that should no longer be needed.
% To do:  
    % DONE: change make_iss.  It will use the f_is_odd value in the list of subIss
    % indices. use the same inv.txt but check for negative colors.
    % bif.txt no longer has negative signs for pitchforks.
    % Add the f_is_odd functionality to isPitchfork.
    %  
    % DONE Put the parameters at the top, and print out the parameters to the
    % command window.
    % Check if I can have f(0) != 0, and input a solution and
    % follow the branches.  For example the constant ISS and us = [1;-4] 
    % (say) when we know that f(1,-4) = 0. I wrote in the nss8.txt
    % manuscript that we can do this, but I haven't tried it.  Line 121 has
    % the starting us vector hard coded

% 
% see sg3_s15 FOR DECIMATION
% Search "tol" for tolerances (small positive floating point numbers)
% Search "cut" or "max" for cutoffs (integers)

% 4/23/2019 note.  I made f not odd so the constant branch has a
% saddle-node bifurcation.  It does not work well.

% 5/7 note:  I am saving sierp_S19_5.m as the version just before I start
% following the branch in the invariant subspace.  I put in a detection of
% pitchfork bifurcations and only followed one branch at a pitchfork.  This
% reduced the number of branches considerably and there are 68 branches
% with this saved version that has smin = -100 and sMax = 130. With smin =
% -200 there are the same number of branches and there don't seem to be any
% bifurcations to the left of -100.  This version does not follow all 5 of
% the branches that bifurcate at s = 75.  -JWS

% Version 6: started 5/7/2019.  The computations are done within the
% invariant subspaces.  The information is read in from the file
% colorLists.txt and then processed and stored in the structure iss.
% (invariant subspaces).  Some of the notation is changed from version 5;
% (1) us is the vector (u, s), replacing w. u = us(1, end-1)
% (2) vs is the vector (v, s)
% (3) j is commonly used to for the invariant subspace W_j with iss(j).
% (4) maxLevel is introduced to limit to primary, secondary, etc. branches
% (5) Special care is taken for the trivial solution.  I tried not to use
% the condition j == 1 but rather use isempty(u)
%
% Future changes (as of 5/17/2019)
% (a) I want to store the various possible schematic functions in the
% structure so we can choose what to plot later.
%
% (b) Include containment data in iss,
% iss(j).contains = {k such that W_k contains W_j}
% Use containment data to reduce the number of branches followed.
% Currently, redundant branches with larger invariant subspaces are
% sometimes followed.
% The containment data is not "rectangular" so I don't know how to read it
% in from a file.  I may just hardcode it and sacrifice the ability to
% include more invariant subspaces at a higher livel.
%
% (c) include the information about pitchfork bifurcations in iss.  This
% can be inclued in containment data.  For example, in iss(j).contains,
% use -k instead of k when there is a pitchform bifurcation from W_j to W_k
%
% (d) The history of a branch should probably include the current bcnt
% For example, the current output says
% Starting branch 62 at s = -8.7153 in invariant subspace 12
% history = [ 1 11 ]
% and it should probably say
% Starting branch 62 at s = -8.7153 in invariant subspace 12
% history = [ 1 11 62 ]
%
% (e) The schematic( ) function is changed for different "y" output.  There
% could/should be a flag at the beginning of the program to choose this.
% Better yet, we should decide on a list of schematic functions to compute
% and save them all to plist and use that flag to choose which to plot.
% The BranchManager could then use that same run to plot different "y"
% functions.
%
% (f) Currently the code works for level n = 2, but should work at higher
% levels provided we only follow the first 5 bifpoints, assuming read in
% the same 13 invariant subspaces, after converting to the higher level.
% It might be possible to just replace colorLists.txt to a new file of a
% matrix of size 13 x 3^n instead of the current 13 x 9
%
% (g) If we want to got to higher level and include bifurcations of the
% trivial solution beyond lam = 125, then the choice of which invariant
% subspaces to include becomes problematic.  At level 2 we can get a
% consistent choice of representatives from each group orbit [W_j].  At
% level 3 I don't think we can do that.  Nandor might be able to do some
% magic, including producing the containment data.
%

function branch(us,vs,j)

global bcnt smin smax blist last_newt_it
global f iss wam print_plist ymax speed
N = size(wam, 1);  %This defines N, which is not a global variable

max_points = 100000; % Default is 100000  (was 1000 before 2023-06-16)
% We used to need this because some brances were infinite loops.
% With MI_list we have eliminating the loops.

usm = us;  % The mother point
stop_branch = false;
points = 0;
% while (smin <= us(end)) && (us(end) <= smax && ~stop_branch && ...
%         points < max_points)
while ~stop_branch
    normvs = norm_usj(vs, j); % standard norm in R^(N+1)
    vs = vs/normvs;
    usg = us + vs * speed; 
    ds_old = vs(end); % Delta s - old
    usn = newton(usg, vs, j);
    vs = usn - us; % new tangent vector
    ds_new = vs(end); % Delta s - old
%     my_norm = norm_usj(vs, j);
%     my_norm = norm(vs);  % Old method.  This is norm in R^m, not R^N
%     vs = vs/my_norm;  % new tangent vector
    
    [MIn, ev_neg_n, ev_pos_n] = morse(usn, j); %OLD
    [MI_list_new, eig_neg_list_new, eig_pos_list_new] = make_MI_list(usn, j);
    if(ds_old*ds_new < 0) && (points > 1)
        if MI_list_old(1) ==  MI_list_new(1)
            fprintf('        !Branch is daughter of a ');
            fprintf('pitchfork bifurcation near last point, ');
            fprintf('stop following it!\n');
            stop_branch = true;
            % This handles
            % MI_list_old(1) = 1,
            % MI_list = 0 because the cutoff.  Should only be 0 at one pt.
            % MI_list_new(1) = 1
            % But there is a bug if 2 or more points give MI = 0
        end
    end
    if points == 0 % MI_list is not defined yet
        dMI_list = 0;
    else
        dMI_list = MI_list_new - MI_list;
    end
    analyzeBif = ~(stop_branch || usn(end) > smax || usn(end) < smin);
    % I am not sure what bugs I had.  I put some code in spawn to make sure
    % the bifpoint is in the window.  I should test this. It looks like I
    % might not be using spawn( ) currently (2020-09-01) but I use
    % addToQueue instead.  I am keeping the spawn( ) function around
    % because it handels the complex pairs of eigenvalues and I eventually 
    % need to handle more than one eigenvalue crossinng the imaginary axis.
    
    % don't analyze bug in a few circumstances.  We might miss a
    % bifurcation in the s-window this way.  Perhaps we should compute the
    % s-value of the bifpoint and include bifpoints in the s-window.
%    if analyzeBif && norm(dMI_list) ~= 0 % MI_list changed
    if ~stop_branch && norm(dMI_list) ~= 0 % MI_list changed
        if dMI_list(1) == 0 && dMI_list(end) == 0
            fprintf('MI_list changed, but not 1st and last components.');
            fprintf('  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
            MI_list % Print these two lists to help understand what happened
            MI_list_new
            % If we very rarely get this message, only compute 1st and last
            % components until we get here and then compute them all.  This
            % will speed up the code significantly
        end
        daughter_ks = find_daughter_ks(dMI_list, j); %NEW
        if abs(dMI_list(1)) == 1
            changingEval = max(MI_list(1), MI_list_new(1));
            if dMI_list(1) == -1
                lam = eig_neg_list(1);
                lam_new = eig_pos_list_new(1);
            else
                lam = eig_pos_list(1);
                lam_new = eig_neg_list_new(1);
            end
            %fprintf('\nCall new dMImother1\n');
            [us0, fold] = dMImother1(j, us, usn, lam, lam_new, changingEval);
            if fold
                blist(bcnt).foldpts = [blist(bcnt).foldpts; ...
                    us0(end), schematic(us0, j)];
                fprintf('MI -> %d at fold point s = %.3g, y = %.3g. ', ...
                    MI_list_new(end), us0(end) , schematic(us0,j) );
                if isempty(daughter_ks)
                    fprintf('\n');
                else
                    fprintf(' Also, daughter(s) in subspace W_');
                    fprintf('%d ', daughter_ks(:));
                    fprintf('\n');
                    subset_i_is_daughter = ismember(iss(j).subsets, daughter_ks);
                    for i = 1: length(iss(j).subsets)
                        if subset_i_is_daughter(i)
                            which_trmd_eig = max(MI_list(i)-MI_list(1),...
                                MI_list_new(i) - MI_list_new(1));
                            k = iss(j).subsets(i);
                            birth_with_fold(j, k, us, usn, ...
                                lam, lam_new, which_trmd_eig);
                        end
                    end
                end
            else  % not a fold, and MI in mother subspace changed by 1
                fprintf('        !Branch is daughter of a ');
                fprintf('transcritical bifurcation near last point, ');
                fprintf('stop following it!\n');
                stop_branch = true;
            end
        elseif abs(dMI_list(1)) == 2
            fprintf('MI within mother subspace changes by 2. Hopf?  ');
            fprintf('No bifurcating branches followed.\n');
        elseif dMI_list(1) == 0
            %The simplest case.
            if length(daughter_ks) > 1
                fprintf('MI -> %d between s = %.3g and s = %.3g with daughters in ', ...
                    MI_list_new(end), us(end), usn(end));
                fprintf('W_%d, ', daughter_ks(1:end-1));
                fprintf('and W_%d.\n', daughter_ks(end));
            else
                fprintf('MI -> %d between s = %.3g and s = %.3g with daughter in W_%d.\n', ...
                    MI_list_new(end), us(end), usn(end), daughter_ks(1));
            end
            % I'm not sure we handled change of MI_list(k) > 1
            %daughter_ks
            subset_i_is_daughter = ismember(iss(j).subsets, daughter_ks);
            for i = 1: length(iss(j).subsets)
                if subset_i_is_daughter(i)
                    k = iss(j).subsets(i);
                    changingEval = max(MI_list(i), MI_list_new(i));
                    if MI_list(i) > MI_list_new(i)
                        birth(j, k, us, eig_neg_list(i), ...
                            usn, eig_pos_list_new(i), changingEval);
                    else %MI_list(i) < MI_list_new(i)
                        birth(j, k, us, eig_pos_list(i), ...
                            usn, eig_neg_list_new(i), changingEval);
                    end
                end
            end
        end
    end % MI_list changed
    if points > 0
        MI_list_old = MI_list;
    end
    MI_list = MI_list_new;
    eig_pos_list = eig_pos_list_new;
    eig_neg_list = eig_neg_list_new;
    if ~stop_branch
        MI = MIn; ev_neg = ev_neg_n; ev_pos = ev_pos_n; us = usn;
    end
    if points == 0
        fprintf('MI =  %d initially.\n', MI);
        plist = full([usm(end), schematic(usm, j), MI, 0]);
        % mother point, with MI of first point
        % do a "reality check"  Is this really a solution?  This is needed
        % because the invariant subspace designation might be wrong.
        u = us(1:end-1);
        if isempty(u)
            uFull = zeros(N,1);
        else
            uFull = iss(j).B*u;
        end
        normg = norm(wam*uFull - f(uFull, us(end)) );
        if(normg > 10^-10) % hardwired tol
            fprintf('Warning! Solution in W_j is not a solution in R^N,');
            fprintf(' normg = %g ', normg);
            fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
        end
    end
    if (smin >= us(end) ) || ( us(end) > smax) ...
            || (abs(schematic(us, j)) > ymax ) || points >= max_points
        stop_branch = true;
    end
    if ~stop_branch  % Should I always increment the points and do the plist, and then set stop_branch?
        points = points + 1;
        plist  = [plist; us(end) schematic(us, j) round(MI, 0) last_newt_it];
    end
end % while loop

usg = us;
sfinal = usg(end);
if sfinal < smin
    usg(end) = smin;  
    remake_last_point;
elseif sfinal > smax
    usg(end) = smax;
    remake_last_point;
%elseif smin < sfinal && sfinal < smax && ~stop_branch
elseif points == max_points
    fprintf('        !Branch stopped with max_points = %i points!  ', ...
        max_points);
    fprintf('Increase speed near line 37, \n\tor consider increasing max_points in branch( ) function.\n');
end % Do none of these if sfinal = smax (as in trivial branch)
    function remake_last_point
        vs = [zeros(length(us)-1,1); 1];
        us = newton(usg,vs,j);
        if exist('MI', 'var') % 2020-08-20 I'm not sure why this is needed
            plist = [plist; us(end) schematic(us, j) MI last_newt_it];
        end
    end

fprintf('\nLast point: s = %.5g, y = %.5g.  uj in W_j and u in R^N follow.\n',...
    us(end), schematic(us, j));
fprintf('uj = ');
fprintf('%.5g ', us(1:end-1));
u = us(1:end-1);
if isempty(u);
    uFull = zeros(N,1);
else
    uFull = full(iss(j).B*u);
end
fprintf('\n u = (');
fprintf('%.5g, ', uFull(1:end-1));
fprintf('%.5g)\n', uFull(end));
%
lam = evals(us, j);
if(norm(imag(lam)) > .1^10)
    fprintf('real parts of ');
end
fprintf('eigenvalues are ');
fprintf('%.2f  ', lam(1:end));
fprintf('\n');
if (print_plist) % output plist at end of branch( )
    fprintf('         s  \t \t y\tMI \t newton its\n');
    [i_max, ~] = size(plist); 
    for i = 1: i_max % using length(plist) = max(height, width) gave a bug
        fprintf('%10g\t%10g\t %d\t %d\n', ...
            plist(i,1), plist(i,2), plist(i,3), plist(i,4));
    end
end

if exist('plist', 'var') % 2020-08-20  I'm not sure why this is needed
    blist(bcnt).pts         = plist(:,1:2);
    blist(bcnt).mi          = plist(:,3);
    blist(bcnt).lastnewt    = plist(:,4);
    blist(bcnt).lastpt      = us;
end

end %%%%%%%% of branch %%%%%%%%%%%%%%%%%%%%%%

function us = newton(us,vs,j)
% input us is the (u, s) guess.  The output is the final us
global verbosity iss f fp dfds last_newt_it

tolNewt = 10^-11;  % 10^-15 is too small, 10^-14 usually works at level 2
maxNewtSteps = 10;

if length(us) == 1
    last_newt_it = 0;
    us = real(us); % needed with with non-symmetric wam
    if verbosity > 1
         fprintf('Found trivial solution with s = %f\n', us(end));
    end
    return % This leaves us unchanged
end
wamj = iss(j).wam;

for it = 1:maxNewtSteps
    
    u   = us(1:end-1);
    s   = us(end);
    G  = [-wamj*u + f(u,s); 0];
    Ju = -wamj + spdiags(fp(u,s),0,length(u),length(u));
    %J   = [Ju, dfds(u,s); vs'*iss(j).dims];
    v = vs(1:end-1);
    vsleft = [ v'*iss(j).dims, vs(end)];
    % vsleft is the row vector associated with vs, from the inner product
    % <u, w> = u' A w in R^(N+1).  The "metric" A = iss(j).dims, with a 1 
    % added at the bottom right of A for the s component.  Note that A is 
    % diagonal and positive definite.  For vs and us in R^N, the
    % constraint in R^N is vs'*us(it) = vs'*usguess, and using vsleft gives
    % the equivalent constraint in the actual us and vs coordinates of W_j.  
    % This proper inner product insures that the first point of a new 
    % daughter branch is not computed to be on the mother branch if the
    % critical eigenvector v is orthogonal to the Mother subspace.    
    J   = [Ju, dfds(u,s); vsleft];

    %J   = [Ju, dfds(u,s); vs']; % old: allows more general s dependence
    %J   = [Ju, u; vs']; % old old version: f(s, x) = s x + nonlinear(x)
    chi = J \ G;
    us   = real(us - chi); % sometimes needed with with non-symmetric wam
    
    if max(abs(chi))<tolNewt%  % was norm(chi)<tolNewt
        break;
    end
end
if it == maxNewtSteps
    fprintf('newton( ) did not converge in %d iterations.  ', ...
        maxNewtSteps);
%     fprintf('Final norm(chi) = %g\n', ...
%         norm(chi));
    fprintf('Final max(abs(chi)) = %g\n', ...
        max(abs(chi)));

else
    if verbosity > 1
        fprintf('Newtons method converged in %d iteration(s) ', it);
        fprintf('to solution with s = %f\n', us(end));
    end
end

last_newt_it = it;

end % of newton( ) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% This computes the smallest  MImax eigenvalues for usj in W_j %%%%%%%%%%%%%%%%%%%%%%%%%%%
function lam = evals(usj, j)
global iss fp wam MImax
MItol = 10^-8; % bug when MItol = 0, and the trivial branch hits s = 75
% with MItol = 10^-10 I got a lot of bogus MI changes on the constant branch
% at level 3
if length(usj) == 1  % trivial solution
    N = size(wam,1);
    u = zeros(N,1);
else
    u = iss(j).B*usj(1:end-1);
end
Jus = -wam + spdiags(fp(u, usj(end)), 0, length(u), length(u));
lam = eigs(-Jus,MImax,'smallestreal');
end % of evals( ) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% this finds list of k in daughter subspaces W_k %%%%%%%%%%%%%%%%%%%%
%
function daughter_ks = find_daughter_ks(dMI_list, j)
global iss
subsets = iss(j).subsets;
daughter_ks = [];
for i=1:length(dMI_list)
    if dMI_list(i) - dMI_list(1) ~= 0
        daughter_ks = [daughter_ks, subsets(i)];
    end
end
% Now, for each daughter k, eliminate the non-k eisslements in
% iss(k).subsets. The MI changes for those subspaces, 
% but there is not a blis bifurcation.
original_daughter_ks = daughter_ks;
for i = 1:length(daughter_ks)
    m = original_daughter_ks(i);
    newBad_ks = setdiff(iss(m).subsets, [m]); 
    daughter_ks = setdiff(daughter_ks, newBad_ks);
end
end % of find_daughter_ks

%%%%%%% This computes the MI_list for usj in W_j %%%%%%%%%%%%%%%%%%%%%%%%%
% For all of the supersets of W_k of W_j (including W_j and the full space,
% compute the Morse index.  Also return the MIth and MI+1st eigenvalue
function [MI_list, eig_neg_list, eig_pos_list] = make_MI_list(usj, j)
global iss fp wam MImax
MItol = 10^-8; % 10^-8 default.  bug when MItol = 0, and the trivial branch hits s = 75
% with MItol = 10^-10 I got a lot of bogus MI changes on the constant branch
% at level 3 (presumably for the PBC graph)
N = size(wam,1);
s = usj(end);
if length(usj) == 1  % trivial solution
    u = zeros(N,1);
else
    u = iss(j).B*usj(1:end-1);
end
subsets = iss(j).subsets;
MI_list = [];
eig_neg_list = [];
eig_pos_list = [];
for k_index = 1: length(subsets)
    k = subsets(k_index);
    [sizeLapk, ~] = size(iss(k).wam);
    if sizeLapk == 0  % this is true for j = k = 1 only !
        MI_list = [0];
        eig_neg_list = [-1];
        eig_pos_list = [1];       
    else
        uk = iss(k).Binv*u;
        Jus = -iss(k).wam + spdiags(fp(uk, s), 0, length(uk), length(uk));
        numEvals = min(MImax, length(uk));
        lam = eigs(-Jus,numEvals,'smallestreal');
        MI = sum(real(lam)<-MItol);
        if MI == 0
            eig_neg = -1; % this should not participate in a bifurcation
            eig_pos = lam(1);
        elseif MI == numEvals % all the eigenvalues are negative
            eig_neg = lam(numEvals);
            eig_pos = 1; % this should not participate in a bifurcation
            
        else
            eig_neg = lam(MI);
            eig_pos = lam(MI+1);
        end
        MI_list = [MI_list, MI];
        eig_neg_list = [eig_neg_list, eig_neg];
        eig_pos_list = [eig_pos_list, eig_pos];
    end
end
end % of make_MI_list( ) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% This computes the full MI for usj in W_j %%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MI, evalneg, evalpos] = morse(usj, j)
global iss fp wam MImax
MItol = 10^-8; % bug when MItol = 0, and the trivial branch hits s = 75
% with MItol = 10^-10 I got a lot of bogus MI changes on the constant branch
% at level 3
if length(usj) == 1  % trivial solution
    N = size(wam,1);
    u = zeros(N,1);
else
    u = iss(j).B*usj(1:end-1);
end
Jus = -wam + spdiags(fp(u, usj(end)), 0, length(u), length(u));
lam = eigs(-Jus,MImax,'smallestreal');
MI  = round(sum(real(lam)<-MItol), 0);
MI  = sum(real(lam)<-MItol);
if MI == 0
    evalneg = -1; % this should not participate in a bifurcation
    evalpos = lam(1);
elseif MI == MImax
    evalneg = lam(MI);
    evalpos = evalneg + 1; % this should not participate in a bifurcation
    
else
    evalneg = lam(MI);
    evalpos = lam(MI+1);
end
end % of morse( ) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sch = schematic(us, j)
global iss wam y
uj = us(1:end-1);
N = size(wam,1);
if isempty(uj)
    sch = 0;
    return
end
% norm_u22 = (uj'*iss(j).dims*uj)/N; % N makes bif diagram converge lev -> inf
% norm_u2 = sqrt(norm_u22);
% u_1 = u_1;

% The cases statement uses the y strings, (search "LaTeX" to find y defs)
switch y
    case '\frac 1 N \|u\|_2^2'
        sch = full(uj'*iss(j).dims*uj)/N;
    case 'u_1'
        sch = full(iss(j).B(1,:)*uj);
    case 'u_N'
        sch = full(iss(j).B(N,:)*uj);
    case 'N u_N'
        sch = N*full(iss(j).B(N,:)*uj);
    case '\frac 1 {\sqrt N} \|u\|_2'
        sch = sqrt(full(uj'*iss(j).dims*uj)/N);
    case '\frac 1 N \|u\|_1'
        absuj = abs(uj);
        sch = full(sum(iss(j).dims*abs(uj)))/N;
    case '\displaystyle{\sum_{i=1}^N \sin(i) x_i}'
        sch = 0;
        ufull = full(iss(j).B*uj);
        for i = 1:N
            sch = sch + sin(i)*ufull(i);
        end
    case '\|u\|_\infty'
        sch = max(abs(full(uj)));
    otherwise
        fprintf('y is not defined. See error message.\n');
end
end % of schematic( )

% We could pass the level n and read different colorLists accordingly 
function make_iss
global iss niss wam graphName auts f_is_odd
% This function populates the structure iss (invariant subspaces) with
% all the information needed concerning the invariant subspaces
%
iss = struct;
% read in partitions as list of colors of invariant subspaces
cd (graphName);
colorLists = dlmread('inv.txt');
subsets = dlmread('bif.txt');
cd '..';
%
[niss, N] = size(colorLists);  % N = number of vertices is not global
% I think we need to make a list of the invariant subspaces when f_is_odd
% is false.  I *think* we want to keep the same numbering, but we could
% suppress the numbers.  There are two ways to proceed:  
% 1: only have the good invariant subspaces in the iss(j).subsets field
% 2: keep the list as it is and have an iss(j).invariant field that is true
% or false.  Using method 1, we need to get that list before we start, or
% perhaps start at the end and do iss(niss) first, and iss(1) last.
% While we are at it, we should probably have a boolean f_of_0_is_0 and if
% this is false not include the (0, 0, ..., 0) invariant subspace.
%      
% Make this list of truly invariant subspaces, and intersect this list with
% the subsets when making the field iss(j).subsets
if f_is_odd
    trulyInvariantSubspaces = [1:niss];
else
    trulyInvariantSubspaces = [];
    for j = [1: niss]
        colorList = colorLists(j,:);
        cmin = min(colorList);
        if cmin >= 0  %  Assumes f(0) = 0, this is the place
            trulyInvariantSubspaces = [trulyInvariantSubspaces, j];
        end
    end
end
for j = [1: niss]
    colorList = colorLists(j,:);
    cmax = max(colorList);
    cmin = min(colorList);
    %     if cmin < 0 %% ~f_is_odd
    %         j
    %         issjInvariant = false
    %     else
    %         issjInvariant = true;
    %     end
    B = sparse(N, cmax);   % Basis matrix for the invariant subspace  W_j
    Binv = sparse(cmax, N);  % Penrose Pseudoinverse
    % I plan to change this to write Pinv at the same time: [P Pinv] output
    for color = [1: cmax]
        rowsPos = find(colorList == color);
        rowsNeg = find(colorList == -color);
        denom = length(rowsPos)+length(rowsNeg);
        for i = [1: length(rowsPos)]
            B(rowsPos(i), color) = 1;
            Binv(color, rowsPos(i)) = 1/denom;
        end
        for i = [1: length(rowsNeg)]
            B(rowsNeg(i), color) = -1;
            Binv(color, rowsNeg(i)) = -1/denom;
        end
    end
    iss(j).B = B;  % sparse integer Basis.  W_j = col(B)
    iss(j).dims = B'*B; % Diagonal matrix of dimensions for each component
    iss(j).Binv = Binv;  % sparse Penrose pseudoinverse of B (non-ints)
    iss(j).P = B*Binv;  % sparse projection matrix onto W_j (non-ints)
    iss(j).wam = Binv*wam*B; % rounding was a bug with noninteger wam
    if cmax == 0  % and the previous loop is
        iss(j).P = sparse(N,N);
    end
    if j <= niss
        subs = subsets(j, 1:nnz(subsets(j,:))); % list read in from bif.txt
        %nnz is number of nonzero elemants
        %only include truly invariant subspaces in list of subsets
        iss(j).subsets = intersect(subs, trulyInvariantSubspaces);
        %         iss(j).subsets = subsets(j, ...
        %             1:nnz(subsets(j,:)));  %nnz is number of nonzero elements
    end
    issjAuts = [];
    for k = 1: length(auts(:,1))
        if isequal(colorList, colorList(auts(k,:)) )
            issjAuts = [issjAuts, k];
        end
        if f_is_odd && isequal(-colorList, colorList(auts(k,:)) )
            issjAuts = [issjAuts, -k];
        end
    end
    iss(j).auts = issjAuts;
    % iss(j).lt % Indices of subspaces less than W_j
end
% sort(abs(iss(2).contains))

end %%% of make_iss %%%%%%%%%%%%%%%%%%%%%

function pitchfork = isPitchfork(jMother, jDaughter, vs)
global iss auts
tol = 10^-5;
vDaughter = vs(1:end-1); % The eigenvector in daughter subspace
v = iss(jDaughter).B*vDaughter;  % The eigenvector in R^n
motherAuts = iss(jMother).auts;  % list of the indices of the automorphisms
numAuts = length(motherAuts);
pitchfork = false;
for i = 2: numAuts
    autIndex = motherAuts(i);  % +/1 i == (auts(i), +/- 1) in Aut(G)xZ_2
    aut = auts(abs(autIndex), :); %
    gv = sign(autIndex)*v(aut);
    if norm(v + gv, Inf) < tol
        pitchfork = true;
        break  % stop the for i = 1: numAuts loop
    end
end
% This does not incorporate f_is_odd yet. 
% I'm not sure about that comment.  I think f_is_odd is incorporated in the
% make_iss function so iss(jMother).auts has signed aut in Aut(G) x Z_2.

% This could be greatly speeded up:  First, does autMother act on
% Wdaughter?  If so, find the action on WDaughter, which is permutation
% combined with a possible sign flip of some coordinates. Then see if that
% group action on vDaughter includes -vDaughter.
% The action on WDaugter will possibly have many duplicates:  That is,
% different elements of autMother give the same action on WDaughter. We
% want to eliminate the redundancy before doing the check.  Alternatively,
% we can just do the check in R^n
% 
% norm(x, 1) = 1-Norm
% norm(x, Inf) = sup norm
end

% call this when MI_list(1) changes by +/- 1.
function [us0, fold] = dMImother1(j, us, us_new, lam, lam_new, maxMI)
vs = us_new - us;
% We are looking for solution for us0 where lam = 0;
x = lam/(lam-lam_new); % The secant method step.
% distance from lam to 0 as a fraction of dist from lam to lam_new
usg = us + x*vs;  % The secant method approximation. We could use Ridders
% This mis-identifies a fold in branch 5 with diamond graph with del = 0.09 in a case
% where lam = 0.1 and lam_new = -5.  I will do at least one step of
% ridders method to see of this helps.
us0 = newton(usg, vs, j);

% lam 
s = us(end);
s0 = us0(end);
% lam_new
s_new = us_new(end);
% if (s0-s)*(s_new-s0) < 0
%     fold = true;
% else
%     fold = false;
% end
fold = (s0-s)*(s_new-s0) <= 0;
end % of dMImotheer1

% The functions colorRGB(j) and colorString(j) need to be coded to agree.
% 2023-06-07 make trivial solution (j = 1) and ISS = R^n (j = niss) black
% For 1< j < niss, cycle through the non-black colors using switch
% I currently use 9 non-black colors for use with the diamond graph
%
% I probably want to change these to more mute colors.  The defaults from
% matlab 2014b are possible
% returns the RGB for the color as a function of the j = issj in W_issj
function rgb = colorRGB(j)
global niss;
if j == niss|| j == 1
    rgb = [0,0,0]; % black
else % mod(j, 9) = 2 for the first nontrivial branch (j = 2)
    switch mod(j, 9) % the j in W_j determines color.
        case 0
            rgb = [1, 0.44, 0.2]; % orange
        case 1
            rgb = [0.5, 0.2, .9]; % plum was [0.4940, 0.1840, 0.5560]
        case 2
            rgb = [1, 0.0, 0.00]; % red
        case 3
            rgb = [.9, .7, 0.1]; % gold was [0.9290, 0.6940, 0.1250]
        case 4
            %rgb = [0, 0.4470, 0.7410];
            rgb = [0, 0.2, 0.9]; % blue
        case 5
            rgb = [0.4660, 0.6740, 0.1880]; % green
        case 6
            rgb = [0.3010, 0.7450, 0.9330]; % cyan
        case 7
            rgb = [.8, 0, .8]; % purple
        case 8
            rgb = [0.7, 0.4, 0.0];  % brown
    end
end
end

% returns the string for the color as a function of the j = issj in W_issj
function string = colorString(j)
global niss;
if j == niss|| j == 1
    string = 'black';
else % mod(j, 9) = 2 for the first nontrivial branch (j = 2)
    switch mod(j, 9) % the j in W_j determines color. Make the two agree!
        case 0
            string = 'orange';
        case 1
            string = 'plumb';
        case 2
            string = 'red';
        case 3
            string = 'gold'; %
        case 4
            string = 'blue'; %
        case 5
            string = 'green'; %
        case 6
            string = 'cyan'; %
        case 7
            string = 'purple'; %
        case 8
            string = 'brown';
    end
end
end

% I'm not sure if I use spawn!

% Can I just find the smallest eigenvalue in each of the
% sub-subspaces and hope that there are abs(MI-MIn) of them?
% MI1 > MI2, ev1 > 0 > ev2
function us_bif = spawn(us1, MI1, ev1, us2, MI2, ev2, j)
global iss wam fp jobs bcnt blist maxLevel maxLevelKickedIn smin smax
tol_zero_eval = 10^-2;  % We can make this smaller after Ridders method
N = size(wam, 1);
vs = us2 - us1;
dimE = MI1 - MI2;
% evGoal is needed for nonsymmetric wam, with nonreal eigenvalues.  We want
% to set the real part of the eigenvalue 0.  I expect ev1 and ev2 to be the
% eigenvalues before and after crossing the imaginary axis.  
% eigs(mat, 1, evGoal) finds the eigenvalues closest to evGoal 
evGoal = 1i*imag(ev1-ev2)/2;  % Pure imaginary (often 0) eigenvalue goal.
if (norm(evGoal) > .1^15)
    fprintf('Next bifurcation appears to be a Hopf.\n');
end
ev1 = real(ev1);
ev2 = real(ev2);
% This should work perfectly for trivial branch, eventually modify this to
% use Ridders Method
x = ev1/(ev1-ev2);  % secant method (no longer used)
x = .5;  % first x for Ridders method
usj_bif = us1 + x*vs;
usj_bif = newton(usj_bif, vs, j);
% Find new smallest eigenvalue:  There might be a problem if a small
% eigenvalue is hanging around.
if length(usj_bif) == 1  % trivial solution
    u_bif = zeros(N,1);
else
    u_bif = iss(j).B*usj_bif(1:end-1);
end
s_bif = usj_bif(end);
Jus = -wam + spdiags(fp(u_bif, s_bif), 0, N, N);
ev_mid = real(eigs(Jus, 1, evGoal)); %eval closest to evGoal (0 or pure imag)
x = .5 + .5*ev_mid/sqrt(ev_mid^2 - ev1*ev2);  % One step of ridder's method
usj_bif = newton(us1 + x*vs, vs, j);
if length(usj_bif) == 1
    u_bif_full = zeros(N,1);
else
    u_bif_full = iss(j).B*usj_bif(1:end-1);
end
s_bif = usj_bif(end);
us_bif = usj_bif;  
if(s_bif <= smin || smax <= s_bif)
    return;
end
%
subsets = abs(iss(j).subsets); % These are the possible daughters
if j == 1
    subsets = setdiff(subsets, [1]); % Don't consider the trivial branch (j=1)
end
while ~isempty(subsets)
    k = subsets(1);
    subsets =setdiff(subsets, [k]);
    u_bif_k = iss(k).Binv*u_bif_full;
    Jusj = -(iss(k).B)'*wam*iss(k).B + iss(k).dims* ...
        spdiags(fp(u_bif_k, s_bif), 0, length(u_bif_k), length(u_bif_k));
    % Jusj = B'wam B is a symmetric matrix, unlike wamj = dims^-1 B' wam B.
    % We  we need to solve the eigenvalues of dims^-1 B' wam B + vector are
    % the generalized eigevectors that satisfy
    % (B' wam B + dims vector)v = lam dims v
    %
    nEvals = min(dimE, length(u_bif_k));
    [evecsk, evalsk] = eigs(-inv(iss(k).dims)*Jusj, nEvals, 'smallestabs');
    %[evecsk, evalsk] = eigs(-Jusj, iss(k).dims, nEvals, 'smallestabs');
    %This caused a bug! With C8directed, trivial solution, and k = 7.  
    %eigs(A, B, 1, "smallestabs') with
    %A = [2 2 0 0; 0 2 2 0; 0 0 2 2; -2 0 0 2]
    %B = diag 2 2 2 2 this caused an error!   
    evalList = diag(evalsk); 
    nZero = sum(abs(evalList) < tol_zero_eval);
    if nZero >= 1 % Take the subsets of W_k out of subsets
        subsets = setdiff(subsets, abs(iss(k).subsets));
        if nZero > 1
            fprintf('non-BLIS bifurcation: %d zero eigenvalues in W_%d.\n',...
                nZero, k);
        end
        if j == k && nZero == 1
            fprintf('Fold point. Found 1 zero eval in W_j of mother.\n');
            return  % Don't spawn: critical eigenvector in mother subspace W_j.
            % Probably a fold point, or maybe a vertical tangent
        end
    end
    [~, levp1] = size(blist(bcnt).hist); % levp1 = level + 1
    lev = levp1 - 1; % The level mother, with first branch lev 0
    if nZero == 1 && k > j && lev >= maxLevel
        maxLevelKickedIn = true;  % don't spawn daughters
    end
    if nZero == 1 && k > j && lev < maxLevel % spawn daughters
        % put vector in job queue
        v = evecsk(:, 1);
        % Now choose the sign consistently.  In W_k, the components of the
        % eigenvector are almost always nonzero, so setting v(1) > 0 works.
        % We often plot full u_1, so the second normalization sets full
        % v(end) > 0 so full v(1) > 0 for even eigenvectors and -1 for odd
        if v(1) < 0
            v = -v;
        end
        fullv = iss(k).B*v;
        if fullv(end)< -10^-8
            v = -v;
        end
        vs = [v; 0];
        us = [u_bif_k; s_bif];
        if smin < s_bif && s_bif < smax && abs(schematic(us, k)) < ymax
            job = length(jobs)+1;
            jobs(job).us = us;
            jobs(job).vs = vs;
            jobs(job).issj = k;
            jobs(job).hist = blist(bcnt).hist;
            if ~isPitchfork(j, k, vs)
                jobs(job+1).us = us;
                jobs(job+1).vs = -vs;
                jobs(job+1).issj = k;
                jobs(job+1).hist = blist(bcnt).hist;
            end
        end        
    end
end

end % of spawn

% Find mother bifpoint in W_j.  We do Newton's method in W_j, where the
% Jacobian is nonsingular.  We want to set an eigenvalue in W_k to 0.
% Then put a new job in the job queue using the bifpoint and critical
% eigenvector in W_k.
% 
% We find the bifpoint using one step of the secant method, or  
% for more accuracy we use the modified regula falsi, or Ridders method. 
% The short explanation is that we want to use a method that is guaranteed
% to bracket the zero eigenvalue.  More than one step of secant method does
% not bracket the root.
%
% This following example shows how we need to be careful.
% Assume the eigenvales, for 0 < x < 1 are -1 + 3x and 1/2, and we have
% computed the solution at x = 0 (evals -1, 1/2, MI = 1) and at 
% x = 1 (evals 1/2 and 2, MI = 0).
% Since we know the eigenvalues, we know that the MI changes at x = 1/3,
% but we want to find this numerically by finding the 0 of some objective
% function of x.
% For a while, we used the 'smallestabs' eigenvalue with eigs, which means
% the objective function, f_abs(x), is defined by
% f_abs(x) = 1/2 for x < 1/6, 
%          = -1 + 3x for 1/6 < x < 1/2
%          = 1/2 for x > 1/2
% so f_abs(x) is not even continuous! This is a bug!
% 
% A better objective function is the max(MI)'th eigenvalue, which in this
% case is the 1st eigenvalue.  We compute the smallest max(MI), and then
% take the largest of these as the objective function.
% Using the 1'st egenvalue with 'smallestreal' in eigs gives
% f_sm(x) = -1 + 3x for x < 1/2
%         = 1/2  for x > 1/2
% f_sm is continuous.
%
% However, using the 'smallestereal' leads to eigs not converging in many
% situations! An example is the current code, where the "graph" is SH5.5
%
% But, doing more that one step of secant method has a bug as follows.  
% If x_1 = 0 and x_2 = 1, then x_3 is the x-intercept of the secant line
% joining (0, -1) and (1, 1/2), which is x_3 = 2/3.
% Then f(x_2) = f(x_3) = 1/2 and the next step of secant method has a
% divide by 0 error.  The secant method gives up the bracketing of the
% root, and this is the most extreme example of the problem.
% I don't think we have mentioned Ridders method in any of our papers. 
% should not happen here because only one eighenvale (the maxMIth one) is
% changing sign. 

function  birth(j, k, us, eig, usn, eign, largerMI)
global iss fp verbosity MImax smin smax ymax
tol_zero_eval = 10^-8;  
max_its = 5;
wamk = iss(k).wam;

Nk = size(wamk,1);

vs = usn - us;
% Use Ridders method to find x in [0, 1] such that the largerMI'th eigenvalue 
% of J restricted to W_k, evaluated at newton(us + x*vs, j), is 0.
% Note that newton(us + x*vs, j)is a point on the mother branch in W_j, and
% the Jacobian J restricted to W_j is nonsingular, so Newton's method works
% fine.
% The eigenvalues and eigenvectors are computed in W_k, the daughter space.

its = 0;
eig1 = eig; x1 = 0;
eig2 = eign; x2 = 1;

% Now do secant method.  Maybe 
more_its = true;
while more_its %(its == 0 || abs(eig2) > tol_zero_eval) && its < 9% secant method 
    its = its+1;
    x3 = x1 + (x2-x1)*eig1/(eig1-eig2); 
    if x3 > 1 +.1^12 || x3 < -.1^12 
        fprintf('x3 = %g which is outside of the rage [0,1].  ', x3);
        fprintf('This is probably an error, becase of the Secant method.');
        fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
    end
    % x3 = (x1 + x2)/2; %This would first step of Ridders method
    us_bif = newton(us + x3*vs, vs, j); % The Jacobian in W_j is invertible
    % when the loop stops us_bif is the bifurcation point
    % now find the eigenvalue in W_k
    u_bif = us_bif(1:end-1); % this is in W_j.  
    if isempty(u_bif)
        u_bif_k = zeros(Nk,1);
    else
        u_bif_k = iss(k).Binv*iss(j).B*u_bif;
    end
    s_bif = us_bif(end);
    Jus = wamk - spdiags(fp(u_bif_k, s_bif), 0, Nk, Nk);
    if false  % Choosing false is the default.  true is quick and dirty
        [evec, eig3] = eigs(Jus, 1, 0); % the one eval closest to 0
    else
        numEvals = min(MImax, Nk);
        [smallestEvecs, smallestEvals] = eigs(Jus, numEvals, 'smallestreal');
        eig3 = smallestEvals(largerMI, largerMI);
        evec = smallestEvecs(:, largerMI);
    end
    %diag(smallestEvals) % Here temporarily.
    % old version
    %[evec, eig3] = eigs(Jus, 1, 0); % the one eval closest to 0
    
    % The next calculations are not used, but they show how eigs doesn't
    % always converge with 'smallestreal' with just a few eigenvalues.  
    % We were able to compute the MI because we never computed so few of the
    % smallest evals.  The next lines define variables that are descriptive
    % in the output but they are not used in the code.
%     closestToZero = eigs(Jus, 1, 0)
%     twoEvalsClosestToZero = eigs(Jus, 2, 0)
%     smallsestEval = eigs(Jus, 1, 'smallestreal')   
%     twoSmallsestEvals = eigs(Jus, 2, 'smallestreal')   
%     largerMIEvalsSmallest = eigs(Jus, changingMI, 'smallestreal')
    % We might be able to use try and catch (matlab commands) to handle 
    % the nonconvergence. For example, try again, asking for more evals.

    eig1 = eig2; x1 = x2;
    eig2 = eig3; x2 = x3;
    more_its = (its < max_its && abs(eig3) > tol_zero_eval);
if verbosity > 1
        if its == 1
            fprintf('               Secant Method eigs');
        end
        fprintf(', %g', eig3);
        if ~more_its
            fprintf('.\n');
        end
    end
end
if abs(eig2) > tol_zero_eval
    fprintf('Secant method did not converge in %d iterations.\n', its)
    fprintf('Final eigenvalue = %.3g whereas tol_zero_eval = %.3g\n', ...
    eig2, tol_zero_eval);
end
us_bif_k = [u_bif_k; s_bif]; % This is the bifurcation point
if verbosity > -1
    %     fprintf('        old point is at s = %f, y = %g\n', ...
    %         us(end), schematic(us, j));
    fprintf('        Bifpoint is at s = %f, y = %g\n', ...
        s_bif, schematic(us_bif_k, k));
    %     fprintf('        new point is at s = %f, y = %g\n\n', ...
    %         usn(end), schematic(usn, j));
    
end
pitchfork = isPitchfork(j, k, [evec; 0]);
if smin <= s_bif && s_bif <= smax  && abs(schematic(us_bif,j)) < ymax
addToQueue(j, us_bif, k, us_bif_k, evec, pitchfork);
if verbosity > 0
    if pitchfork
        fprintf('        pitchfork bifurcation to W_%d: evec = (', k);
    else
        fprintf('        transcritical bifurcation to W_%d: evec = (',k);
    end
    %fprintf('        Critical eigenvector in W_%d is (', k);
    if length(evec)> 1
        fprintf('%.4f, ', evec(1:end-1));
    end
    fprintf('%.4f)\n', evec(end));
end
end
%fprintf('Mother in W_%d gives birth to daughter in W_%d\n', j, k);
% Find mother in W_j, using one step of Secant method or Secant method or
% Ridders method.  The peoblem that caused us to go to Ridders method
% should not happen here because only one eighenvale (the maxMIth one) is
% changing sign. 
% Then find critical eigenvector in W_j and put a new job in the queue.
end % birth( )

% mother branch only comes in through j, so we can compute isPitchfork
% us is daughter in W_k x R, and evec is in W_k  
% Add a job to the jobs queue if the level does not exceed maxLevel
function addToQueue(j, usj, k, usk, evec, pitchfork)
global blist bcnt maxLevel maxLevelKickedIn iss jobs
% bcnt is the number of the branch that is currently being followed
blist(bcnt).bifpts = [blist(bcnt).bifpts; usj(end),schematic(usj,j)];
[~, levp1] = size(blist(bcnt).hist); % levp1 = level + 1
lev = levp1 - 1; % The level mother, with first branch lev 0
if lev >= maxLevel
    maxLevelKickedIn = true;  % don't add to queue. Flag used in main
else
    % put vector in job queue
    % Now choose the sign consistently.  In W_k, the components of the
    % eigenvector are almost always nonzero, so setting v(1) > 0 works.
    % We often plot full u_1, so the second normalization sets full
    % v(end) > 0 so full v(1) > 0 for even eigenvectors and -1 for odd
    if evec(1) < 0
        evec = -evec;
    end
    fullv = iss(k).B*evec;
    if fullv(end)< -10^-8
        evec = -evec;
    end
    vs = [evec; 0];
    job = length(jobs)+1;
    jobs(job).us = usk;
    jobs(job).vs = vs;
    jobs(job).issj = k;
    jobs(job).hist = blist(bcnt).hist;
    %    if ~isPitchfork(j, k, vs)
    if ~pitchfork
        jobs(job+1).us = usk;
        jobs(job+1).vs = -vs;
        jobs(job+1).issj = k;
        jobs(job+1).hist = blist(bcnt).hist;
    end
end
end % addToQueue

% This function is called whith a fold with MI_list(1) changing by 1 (note
% that the index of the mother branch is always 1)
% and MI_list(index of k) - MI_list(1) changing by 1
% This differs from birth (without fold) in that here the eigenvalue is in
% W_j (the mother) and there it was in W_k (the daughter)
function birth_with_fold(j, k, us_o, us_n, eig_jo, eig_jn, which_trmd_eig)
global iss fp verbosity MImax smin smax ymax

tol_zero_eval = 10^-8;
max_secant_its = 7;

wamj = iss(j).wam;
wamk = iss(k).wam;

Nj = size(wamj, 1);
Nk = size(wamk, 1);

vs = us_n - us_o;
% we use secant method to find x in [0, 1] such that newton(us + x*vs, j)
% has a zero eigenvalue.  Since a single eigenvalue is near zero, the
% secant method should work OK, unlike what forced us to go to ridders
% method when we had several eigenvalues.

%its = 0;
x1 = 0; [eig1, ~, ~, ~] = eig_k_of(0);
x2 = 1; [eig2, ~, ~, ~]= eig_k_of(1);

secant_its = 0;
more_its = true;
while more_its
    secant_its = secant_its+1;
    x3 = x1 + (x2-x1)*eig1/(eig1-eig2); % One step of secant method
    if x3 > 1 || x3 < 0 
        fprintf('x3 = %g which is outside of the rage [0,1].  ', x3);
        fprintf('This might be, becase of the Secant method.');
        fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
    end
    [eig3, us_j, us_k, J_k] = eig_k_of(x3);
    more_its = (secant_its < max_secant_its) && abs(eig3) > tol_zero_eval;
    if verbosity > 1
        if secant_its == 1
            fprintf('        Secant Method eigs');
        end
        fprintf(', %g', eig3);
        if ~more_its
            fprintf('.\n');
        end
    end
    eig1 = eig2; eig2 = eig3;
    x1 = x2; x2 = x3;
end
if abs(eig2) > tol_zero_eval
    fprintf('Secant method did not converge in %d iterations.\n', secant_its)
    fprintf('Final eigenvalue = %.3g whereas tol_zero_eval = %.3g\n', ...
    eig2, tol_zero_eval);
end

% Now the bifpoint us_k, the near-zero eigenvector eig3,
% and J_k, the Jacobian in W_k are stored. compute the critical eigenvector
% using 'smallestabs' in eigs)( ).  The eigenvalue should be eig3.  If not,
% "Something is wrong!!!" This is
% perhaps inefficient to compute the eigenvector again, but the way we
% eliminate the eigenvalues of J_k that were also eigenvalues of J_j makes
% this the easiest way to find the critical eigenvector.

[evec_k, eig_k] = eigs(J_k, 1, 'smallestabs');
if abs(eig_k - eig3) > .1^10
    fprintf('Something is wrong.!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
end
% This next bit is needed because sometimes the are multiple eigenvalues
% with this same eigenvalue and want the eigenvector orthogonal to W_j
proj_evec_k_onto_j = iss(k).Binv*iss(j).P*iss(k).B*evec_k;
if norm(proj_evec_k_onto_j) > .1
    %fprintf('orthogonal projection');
    evec_k = evec_k - proj_evec_k_onto_j;
    %fprintf('normalized orthogonal projection');
    if norm(evec_k) < .1^5 % This shouldn't happen
        fprintf('Danger, Will Robinson!!!!!!!!!!!!!!!!\n');
        evec_k = proj_evec_k_onto_j;
    else
        evec_k = evec_k/norm(evec_k);
    end
end

if verbosity > -1
%     fprintf('               old point is at s = %f, y = %g\n', ...
%         us_o(end), schematic(us_o, j));
    fprintf('        Bifpoint is at s = %f, y = %g\n', ...
        us_k(end), schematic(us_k, k));
%     fprintf('               new point is at s = %f, y = %g\n\n', ...
%         us_n(end), schematic(us_n, j));
end

pitchfork = isPitchfork(j, k, [evec_k; 0]);
if smin <= us_k(end) && us_k(end) <= smax && abs(schematic(us_k, k)) < ymax
    addToQueue(j, us_j, k, us_k, evec_k, pitchfork);
    if verbosity > 0
        if pitchfork
            fprintf('        pitchfork bifurcation to W_%d: evec = (', k);
        else
            fprintf('        transcritical bifurcation to W_%d: evec = (',k);
        end
        if length(evec_k)> 1
            fprintf('%.4f, ', evec_k(1:end-1));
        end
        fprintf('%.4f)\n', evec_k(end));
    end
end
    function us_of_x = us_of(x)
        if x == 0
            us_of_x = us_o;
        elseif x == 1
            us_of_x = us_n;
        else
            us_g = us_o + x* vs;
            us_of_x = newton(us_g, vs, j);
        end
    end

    function [eig, us_j, us_k, J_k] = eig_k_of(x)
        us_of_x = us_of(x);
        u_of_x = us_of_x(1:end-1); % this is in W_j.
        if isempty(us_of_x)
            u_k_of_x = zeros(Nk,1);
        else
            u_k_of_x = iss(k).Binv*iss(j).B*u_of_x;
        end
        s_of_x = us_of_x(end);
        us_k = full([u_k_of_x; s_of_x]);
        J_j = wamj - spdiags(fp(u_of_x, s_of_x), 0, Nj, Nj);
        eig_j = eigs(J_j, 1, 0); % the one eval closest to 0
        numEvals = min(MImax, Nj);
        lam_jSmallest = eigs(J_j, numEvals,'smallestreal');
        
        J_k = wamk - spdiags(fp(u_k_of_x, s_of_x), 0, Nk, Nk);
        numEvals = min(MImax, Nk);
        lam_kSmallest = eigs(J_k, numEvals, 'smallestreal');
        %x
%         eigs_k = diag(eigs_k);  %Now eigs_k is a list, before it was matrix
        trimmed_lam_k = lam_kSmallest;
        iMax = sum(lam_jSmallest < lam_kSmallest(end) + .1^5);
        for i = 1: iMax
            beforeLength = length(trimmed_lam_k);
            trimmed_lam_k = trimmed_lam_k(abs(trimmed_lam_k-lam_jSmallest(i))>.1^12);
            afterLength = length(trimmed_lam_k);
            if beforeLength - afterLength == 2
                trimmed_lam_k = sort([trimmed_lam_k; lam_jSmallest(i)]);
            end
        end
        
%        b = a(a~=3); % b is the vector a with the elements 3 removed
%         tol = .1^10; % tolerance for eigenvalues so close 
%                      % that eigenvectors get mixed up
%         if abs(eigs_k(2) - eig_j ) > abs(eigs_k(1) - eig_j) + tol
%             eig = eigs_k(2);
%             crit_evec = evecs_k(:,2);
%         elseif abs(eigs_k(2) - eig_j ) < abs(eigs_k(1) - eig_j) - tol
%             eig = eigs_k(1,1);
%             crit_evec = evecs_k(:,1);
%         else % eigenvalues are so close that critical eigenspace is 2-D
%             eig = eigs_k(1);
%             proj_evec1 = iss(k).Binv*iss(j).P*iss(k).B*evecs_k(:,1);
%             proj_evec2 = iss(k).Binv*iss(j).P*iss(k).B*evecs_k(:,2);  
%             if norm(proj_evec1) < norm(proj_evec2)
%                 crit_evec = evecs_k(:,1) - proj_evec1;
%             else
%                 crit_evec = evecs_k(:,2) - proj_evec2;
%             end
%             crit_evec = (1/norm(crit_evec))*crit_evec;
%         end
        us_j = us_of_x;
        %eig = trimmed_lam_k(which_trmd_eig,which_trmd_eig);
        eig = trimmed_lam_k(which_trmd_eig);
    end % function eig_k_of(x)
end % function birth_with_fold

% The next function controls how long the vs tangent vector is in the
% tGnga.  The tangent vector always has my_norm = 1, which means that the
% regular norm in R^(N+1) = speed.  As you can see from this bit of code I
% experimented with various norms in R^N and also with having different
% speed in the s direction and the u direction.  I decided that one speed
% is OK with a variety of systems, after I changed from switch 0 (the norm
% in R^m = W_m, the mother subspace, to switch 1 (the norm in R^N).  I left
% the code in this complicated form to facilitate re-introduction of two
% different speeds.  To go back to this, define speed_s and speed_u in the
% main function, and make them global there and here.  The output of speed
% should be modified too.
%
% speed_s and speed_u replace del on 2020-09-09
% Default is speed_u = speed_s
% speed replaces speed_u = speed_s on 2020-09-19
% norm of u in R^N replaces regular norm of u on 2020-09-19
% 
function my_norm = norm_usj(us, j)
global iss
u = us(1:end-1); 
s = us(end);
if isempty(u)
    norm_u2 = 0;
else
    N = sum(full(diag(iss(j).dims))); % I should maybe make N global
    norm_u2 = u'*iss(j).dims*u; % this reproduces norm in R^N
%     switch 1 % I think switch 1 is good for SH5.5 (N=60) and for SG9 (N=9, 
%              % with scaling in the direction, 
%              % and for diamond (N=4) with speed_s = speed_u 
%              % So I am planning to use just one speed and the norm in R^n
%         case 0 % What we did until recently
%             norm_u2 = norm(us)^2;  % in R^m (m = dim(W_m)), not in R^N
%         case 1
%             norm_u2 = u'*iss(j).dims*u; % this reproduces norm in R^N
%         case 2
%             norm_u2 = u'*iss(j).dims*u/N;
%         case 3
%             norm_u2 = (sum(iss(j).dims*abs(u))^2)/N;%||u||_1 taxicab in R^N
%         case 4
%             norm_u2 = max(u.^2);
%     end
end
my_norm = sqrt(norm_u2 + s^2); % This is the standard norm in R^(N+1)
% speed_s = speed;
% speed_u = speed;
% my_norm = sqrt(norm_u2/speed_u^2 + s^2/speed_s^2); % This is the new norm
end % function norm_usj