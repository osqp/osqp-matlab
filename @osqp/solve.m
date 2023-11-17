%%
function varargout = solve(this, varargin)
    % SOLVE solve the QP

    nargoutchk(0,1);  %either return nothing (but still solve), or a single output structure
    [out.x, out.y, out.prim_inf_cert, out.dual_inf_cert, out.info] = osqp_mex('solve', this.objectHandle);
    if(nargout)
        varargout{1} = out;
    end
    return;
end