function osqp_sfun(block)
% Level-2 MATLAB file S-Function for OSQP.

  setup(block);
  
%endfunction

function setup(block)
  
  % basic problem sizes
  m  = block.DialogPrm(1).Data;
  n  = block.DialogPrm(2).Data;

  %% Register number of input and output ports
  block.NumInputPorts  = 9;  %refactor,warmstart,Px,Px_idx,Ax,Ax_idx,q,l,u
  block.NumOutputPorts = 17; %x,y,prim_inf_cert,dual_inf_cert, various info values
  
  %outputPortSize = [n,m,n,m,ones(1,14)];
  
  % Register the parameters.
  block.NumDialogPrms     = 9; %m,n,P,A,q,l,u, emptysolver, options
  block.DialogPrmsTunable = repmat({'Tunable'},[1 block.NumDialogPrms]);

  %% Setup functional port properties to dynamically
  %% inherited.
  block.SetPreCompInpPortInfoToDynamic;
  block.SetPreCompOutPortInfoToDynamic;
 
  %Set input port properties
  for i = 1:block.NumInputPorts
    block.InputPort(i).DirectFeedthrough = true;
    block.InputPort(i).DatatypeID  = 0;  % double
    block.InputPort(i).Complexity  = 'Real';
  end
  
  %Set output port properties

  %initialise all sizes to scalar
  outdims = ones(1,block.NumOutputPorts);
    
  %The first four output ports are x,y,prim_cert, dual_cert.   
  %everything else is scalar
  outdims(1:4) = [n m m n];
    
  for i = 1:block.NumOutputPorts
    block.OutputPort(i).DatatypeID  = 0;  % double
    block.OutputPort(i).Complexity  = 'Real';
    block.OutputPort(i).SamplingMode = 'Sample';
    block.OutputPort(i).Dimensions = outdims(i);
  end
  
  %% Set block sample time to inherited
  block.SampleTimes = [-1 0];
  
  %% Set the block simStateCompliance to default (i.e., same as a built-in block)
  block.SimStateCompliance = 'DefaultSimState';

  %% Run accelerator on TLC
  block.SetAccelRunOnTLC(true);
  
  %% Register methods
  block.RegBlockMethod('Outputs',@Output);  
  block.RegBlockMethod('PostPropagationSetup', @DoPostPropSetup);
  block.RegBlockMethod('SetInputPortSamplingMode', @SetInpPortFrameData);
  block.RegBlockMethod('SetInputPortDimensions', @SetInpPortDims);
  
%endfunction

function SetInpPortDims(block, idx, di)
  block.InputPort(idx).Dimensions = di;
%endfunction

function DoPostPropSetup(block)

  %solver setup
  
  %the solver object is created as a parameter,
  %since there is seemingly no place to stash it in the 
  %internal memory of the block.   In C it would be possible
  %to store it in pwork, but there's no .m s function analogy
%   m       = block.DialogPrm(1).Data;
%   n       = block.DialogPrm(2).Data;
  P       = block.DialogPrm(3).Data;
  A       = block.DialogPrm(4).Data;
  q       = block.DialogPrm(5).Data;
  l       = block.DialogPrm(6).Data;
  u       = block.DialogPrm(7).Data;
  solver  = block.DialogPrm(8).Data;
  opts    = block.DialogPrm(9).Data;
  
  if(iscell(opts))
    solver.setup(P,q,A,l,u,opts{:})
  else
    solver.setup(P,q,A,l,u,opts)
  end

%endfunction



function Output(block)

  %grab all of the external signals
  %-----------------------------------------------------
  refactor  = block.InputPort(1).Data;
  warmstart = block.InputPort(2).Data;
  Px        = block.InputPort(3).Data;
  Px_idx    = block.InputPort(4).Data;
  Ax        = block.InputPort(5).Data;
  Ax_idx    = block.InputPort(6).Data;
  q         = block.InputPort(7).Data;
  l         = block.InputPort(8).Data;
  u         = block.InputPort(9).Data;
  
  
  %grab parameters and solver data
  %-----------------------------------------------------
  m      = block.DialogPrm(1).Data;
  n      = block.DialogPrm(2).Data;
  solver = block.DialogPrm(8).Data;
  
  
  %construct data updates lists
  %-----------------------------------------------------
  %a list of updates to push to the solver
  updates = {};
  
  %if refactoring is enabled, then add updates to Px etc
  if(refactor)
      updates = [updates,{'Px',Px,'Px_idx',Px_idx,'Ax',Ax,'Ax_idx',Ax_idx}];
  end
  
  %add updates to q/l/u if they are not NaN valued
  if(~any(isnan(q))), updates = [updates,{'q',q}]; end
  if(~any(isnan(l))), updates = [updates,{'l',l}]; end
  if(~any(isnan(u))), updates = [updates,{'u',u}]; end
 
  if(length(updates) > 0)
      solver.update(updates{:})
  end
  
  %reset variables to zero on cold start
  %-----------------------------------------------------
  if(~warmstart)
      solver.warm_start('x',zeros(n,1),'y', zeros(m,1));
  end
  
  %solve and map to outputs
  %-----------------------------------------------------
  sol = solver.solve();
  %will be gathered in first bus
  block.OutputPort(01).Data = sol.x;
  block.OutputPort(02).Data = sol.y;
  block.OutputPort(03).Data = sol.prim_inf_cert;
  block.OutputPort(04).Data = sol.dual_inf_cert;
  
  %will be gathered in second bus
  block.OutputPort(05).Data = sol.info.iter;  
  block.OutputPort(06).Data = sol.info.status_val;
  block.OutputPort(07).Data = sol.info.status_polish;
  block.OutputPort(08).Data = sol.info.obj_val;
  block.OutputPort(09).Data = sol.info.pri_res;
  block.OutputPort(10).Data = sol.info.dua_res;
  block.OutputPort(11).Data = sol.info.setup_time;
  block.OutputPort(12).Data = sol.info.solve_time;
  block.OutputPort(13).Data = sol.info.update_time;
  block.OutputPort(14).Data = sol.info.polish_time;
  block.OutputPort(15).Data = sol.info.run_time;
  block.OutputPort(16).Data = sol.info.rho_updates;
  block.OutputPort(17).Data = sol.info.rho_estimate;
  
%endfunction

function SetInpPortFrameData(block, idx, fd)
  
  block.InputPort(idx).SamplingMode = fd;
  block.OutputPort(1).SamplingMode  = fd;
  
%endfunction
