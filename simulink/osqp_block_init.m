function osqp_block_init(parentBlock)

    %------------------------------------------------------
    % Extract DialogParameters from the block
    refactor_style              = get_param(parentBlock,'refactor_style');
    warmstart_style             = get_param(parentBlock,'warmstart_style');
    has_external_q_input        = get_param(parentBlock,'has_external_q_input');
    has_external_bound_input    = get_param(parentBlock,'has_external_bound_input');


    %------------------------------------------------------
    %Handle the configurable input ports that use pulldowns
    %------------------------------------------------------

    thePortnames         = {'refactor_ext','warmstart_ext'};
    thePulldownValues    = {refactor_style, warmstart_style};

    for i = 1:length(thePortnames)

        portname = thePortnames{i};
        pdVal    = thePulldownValues{i};

        %Configure the block with name = <portname> to be
        %of either constant or inport type
        %------------------------------------
        switch pdVal

            case {'Never','Always'}
                osqp_replace_block(parentBlock,portname,'Constant')

            case {'Triggered'}        %From external signal
                osqp_replace_block(parentBlock,portname,'Inport')

            otherwise
                error('Unrecognized pulldown value');

        end % switch

        %For constant cases, set an appropriate disable/enable signal
        %------------------------------------
        switch pdVal

            case 'Never'
                set_param([parentBlock, '/' portname],'Value','0')

            case 'Always'
                set_param([parentBlock, '/' portname],'Value','1')

            case 'Triggered'
                %do nothing, value is from external port

            otherwise
                error('Unrecognized pulldown value');

        end %switch

    end %end loop

    %------------------------------------------------------
    %Enable/Disable external matrix data if refactors are enabled
    %------------------------------------------------------

    thePortnames         = {'P_idx_ext','P_val_ext','A_idx_ext','A_val_ext'};

    %Configure the block with name = <portname> to be
    %of either constant or inport type
    %------------------------------------
    for j = 1:length(thePortnames)
        portname = thePortnames{j};

        switch refactor_style

            case {'Always','Triggered'}
                osqp_replace_block(parentBlock,portname,'Inport')

            case 'Never'
                osqp_replace_block(parentBlock,portname,'Constant')
                set_param([parentBlock, '/' portname],'Value','NaN')

            otherwise
                error('Unrecognized pulldown value');

        end


    end %end loop


    %------------------------------------------------------
    %Handle the configurable input ports that use tickboxes
    %------------------------------------------------------

    thePortGroups    = {{'q_ext'},{'l_ext','u_ext'}};
    theTickValues    = {has_external_q_input,has_external_bound_input};

    for i = 1:length(theTickValues)

        portnames = thePortGroups{i};
        tickVal    = theTickValues{i};

        %Configure the block with name = <portname> to be
        %of either constant or inport type
        %------------------------------------
        for j = 1:length(portnames)
            portname = portnames{j};

            switch tickVal

                case 'off'   %unticked
                    osqp_replace_block(parentBlock,portname,'Constant')
                    set_param([parentBlock, '/' portname],'Value','NaN')

                case 'on'   %ticked
                    osqp_replace_block(parentBlock,portname,'Inport')

                otherwise
                    error('Unrecognized dialog value');

            end %switch

        end % for j

    end %end loop

end %end function








function osqp_replace_block(parentBlock,portname,newtype)

    %replace the named block with one of the specified type, but
    %only do this if the new type is different from the current one

    %this is necessary since replacing an input port block with
    %the same type has the effect of disconnecting the port
    
    oldtype = get_param([parentBlock, '/' portname],'BlockType');

    if(~strcmp(oldtype,newtype))
        replace_block(parentBlock,'Name',portname,newtype,'noprompt');
    end


end %end function
