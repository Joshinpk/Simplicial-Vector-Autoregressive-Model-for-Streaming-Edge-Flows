classdef TirsoState
    % TirsoState objects encode the state of a TIRSO online estimator 
    %            this class will be used together with the Tirso class.
    
properties
    m_buffer % P  x N  matrix % contains the P previous vectors y[t-p]
                              % in reverse order (newest up)
    m_A      % PN x N  matrix
    m_Phi    % PN x PN matrix
    m_R      % PN x N  matrix % collects the N PN-vectors r_n
    t        % keeps track of the time instant (for diminishing stepsize)
    % Total memory: PN*(P+2)*N + P*N positions
    % = (PN*(P+2)+1)*N positions, O(P^2*N^2)
end

methods
    function obj = TirsoState(ch_message)
        % constructor. Do not use directly. Construct a Tirso and call its
        % initialize method instead.
        errormsg = ['The correct way of constructing a TirsoStateobject'... 
                    ' is calling the initialize method of a Tirso object'];
        if not(exist('ch_message', 'var'))
            error(errormsg)       
        else
            switch ch_message
            case 'calling from Tirso.initialize'
                % OK
            otherwise
                warning (errormsg);
            end
        end
    end
    function ts_out = updateBuffer(obj, v_y)
        % updateBuffer: shift the buffer one position and inserts 
        % a new data vector.
        % Input:  v_y    N-vector 
        % Output: ts_out object encoding new state 
        %
        % Note: TirsoState is not handle; in other words:
        % updateBuffer does not modify the object over which is called
        assert( iscolumn(v_y) && length(v_y)==size(obj.m_buffer,2) );
        ts_out = obj;
        ts_out.m_buffer(2:end, :) = obj.m_buffer(1:end-1, :);
        ts_out.m_buffer(1, :) = v_y;
    end
    
    function [b_out, message] = isConsistentWith(obj, tirso_in)
        % isConsistentWith: check whether this object is consistent with
        % its corresponding instance of Tirso.
        % Input:  tirso_in Object of class Tirso
        % Output: b_out   boolean valued 1 when consistent (OK), 0 when not
        %         message char    telling what was wrong
        assert(isa(tirso_in, 'Tirso'));
        P = tirso_in.order;
        N = tirso_in.noOfNodes;
        
        b_out = false;
        if not(    isequal(size(obj.m_buffer), [P, N]))
            message = 'Incorrect size of buffer';
        elseif not(isequal(size(obj.m_A),      [P*N, N]))
            message = 'Incorrect size of A matrix';
        elseif not(isequal(size(obj.m_R),      [P*N, N]))
            message = 'Incorrect size of R matrix';
        elseif not(isequal(size(obj.m_Phi),    [P*N, P*N]))
            message = 'Incorrect size of Phi matrix';
        else
            b_out = true; message = '';
        end
    end
    
    function t3_out = t_estimatedCoefficients(obj)
        % t_estimatedCoefficients: express estimated coefficients as a 3-way
        % tensor. The VAR parameters in the m_A property are encoded 
        % Output: t3_out N x N x P tensor encoding VAR parameters
        
        [P, N] = size(obj.m_buffer);
        t3_out = permute(reshape(obj.m_A, [P, N, N]), [3 2 1]); 
    end
    
    function m_out = m_estimatedAdjacency(obj)
        % m_estimatedAdjacency express: weighted adjacency matrix from the 
        % estimated VAR parameters. The weights equal the 2-norm of the
        % P-vector containing the impulse response corresponding to each
        % edge.
        % Output:m_out N x N matrix weighted adjacency
        
        [P, N] = size(obj.m_buffer);
        m_A_fat = reshape(obj.m_A, [P N*N]);
        v_norms = sqrt(sum(m_A_fat.^2));
        m_out = reshape(v_norms, [N N])';
    end
        
    function v_out = predictFromBuffer(obj)
        % predictFromBuffer: performs 1-step ahead prediction, just by
        % generating the v_g vector from the buffer and multiplying m_A*v_g
        % Output: N-vector with predicted values
        
        v_g = obj.m_buffer(:);
        v_out = obj.m_A'*v_g;
    end
    
    function m_out = predictMany(obj, m_pastValues, K_toPredict)
        % predictMany predict an arbitrary (k) number of steps ahead
        % Input: m_pastValues P x N matrix containing P past value vectors  
        %                                  (oldest up)
        %        K_toPredict  scalar       number k of steps to predict  
        % Output: m_out       K x N matrix containing predictions for all
        %                     time horizons
        [P, N] = size(m_pastValues);        
        assert(isequal(size(obj.m_A), [P*N, N]));
        
        ts_local = obj;
        ts_local.m_buffer = flipud(m_pastValues); 
        % values stored in buffer in reverse time order
        
        m_out = zeros(K_toPredict, N);
        for k = 1:K_toPredict
            m_out(k, :) = ts_local.predictFromBuffer;
            ts_local = ts_local.updateBuffer(m_out(k,:)');
        end
    end
    
    function m_out = predictManyFromBuffer(obj, K_toPredict)
        % predictManyFromBuffer predicts k steps ahead taking the current
        % contents of the buffer as the past data
        % Input: K_toPredict  scalar       number k of steps to predict
        
        m_out = obj.predictMany(flipud(obj.m_buffer), K_toPredict);
    end
    
end

methods (Static)
    function m_A = get_m_A(t3_in)
        % get_m_A translate set of P square matrices with VAR parameters
        % into fat matrix (such that m_A * v_g is the prediction).
        % Input:  t3_in N x N x P tensor containing P square matrices
        % Output: m_A   P*N x N   matrix such as the one stored in the 
        %                         m_a property
        %
        % Note: this is the inverse operation of the method
        % t_estimatedCoefficients
        
%         assert(ndims(t3_in) == 3);
        [N_out, N_in, P] = size(t3_in);
        assert(N_out==N_in); 
        N = N_in;
        m_A = reshape(permute(t3_in, [3 2 1]), [P*N, N]);
    end
end
    
end