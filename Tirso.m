classdef Tirso
    % TIRSO Efficient implementation of [zaman2019online]
properties
    regPar           % lambda (can be a scalar or a vector)
    forgettingFactor % gamma
    order            % P
    noOfNodes        % N
    h_stepsize       = @(ts)0.5/eigs(ts.m_Phi, 1);
    
    b_shrinkSelfLoops = 0; % 0:Bolstad; 1:GroupLasso

    % possible improvement: add intercept
end

methods
    function ts = initialize(obj, sigma2_init, m_first_samples, t_A)
        % initialize produces a TirsoState object with the initial values
        % Input: sigma2_init initial value for the diagonal of m_Phi
        %        m_first_samples P x N matrix containing P first samples
        %                                     (oldest up)
        %        t_A (optional)  N x N x P (3-way) tensor containing
        %                                     initial VAR parameters
        % Output: ts             object of class TirsoState
        ts = TirsoState ('calling from Tirso.initialize');
        P = obj.order;
        N = obj.noOfNodes;
        assert(isequal(size(m_first_samples), [P, N]), ...
            'size of m_first_samples should be P x N (obj.order -by- obj.noOfNodes')
        if isscalar(sigma2_init)
            ts.m_Phi = sigma2_init*eye(P*N);
        elseif ismatrix(sigma2_init)
            ts.m_Phi = sigma2_init;
        else
            error 'wrong format sigma2'
        end
        ts.m_R   = zeros(P*N,N); % r_n ^= ts.m_R(:,n)
        
        %ts.m_A = zeros(P*N, N);
        if exist('t_A', 'var')
            assert(isequal(size(t_A), [N, N, P]), ...
                'inconsistent size of initial parameter tensor');
        else
            t_A = zeros(N, N, P);
            t_A(:,:, 1) = 0.99*eye(N);
        end
        ts.m_A = ts.get_m_A(t_A); % turn it into a P*N x N matrix
        
        ts.m_buffer = flipud(m_first_samples);
        % we flipud because the input matrix is sorted from early to late,
        % but we store the samples in the buffer in reverse order.
        
        %ts.maxEigenvalue_Phi = sigma2_init;
        assert(ts.isConsistentWith(obj));
        % TODO: maybe also keep an (NxN) sparsity map of A
    end

    function [ts_out, alpha, m_B] = update(obj, ts_in, v_y)
        % update Performs the Tirso estimator update
        % Input:  ts_in object of class TirsoState (remember to initialize)
        %         v_y   N- vector containing the current data vector
        % Output: ts_out object of class TirsoState encoding the new state
        %         alpha  scalar   the step size used for the update
        %         m_B    auxiliary matrix for advanced use. Can be ignored.
        %
        % Note: TirsoState is not handle, so this call does not modify
        % ts_in.
        % Example: ts = my_tirso.update(ts, m_y(:,t))
        [b_c, message] = ts_in.isConsistentWith(obj);
        assert(b_c, message)
        
        % Receive data vector y[t]
        assert(isequal(size(v_y), [obj.noOfNodes, 1]));
        
        % Form g[t] via (4)
        % % g[t] = vec([y[t-1], ..., y[t-P]]') (4)
        v_g = ts_in.m_buffer(:); % g is the vec of a P*N matrix
        
        % m_Phi[t] = gamma*m_Phi[t-1] + mu*g[t]*g[t]'
        m_Phi = obj.forgettingFactor  * ts_in.m_Phi + ...
             (1-obj.forgettingFactor) * (v_g*v_g');
         
        % With these lines I tried to make it faster but did not gain
        % anything!
%         m_summand1 = obj.forgettingFactor  * ts_in.m_Phi;
%         v_factor1 = sqrt(1-obj.forgettingFactor)*v_g;
%         m_summand2 = v_factor1*v_factor1';
%         m_Phi = m_summand1 + m_summand2;
     
        
        % for n = 1:N
            % r_n[t] = gamma * r_n[t-1] + mu * y_n[t] * g[t]
            % v_n[t] = m_Phi[t] * a_n[t] - r_n[t];
            % These can be written in matrix form! we take this first
            % assignments out of the for loop:
        % end
        m_R =  obj.forgettingFactor  * ts_in.m_R + ...
            (1-obj.forgettingFactor) * v_g*v_y';
        m_V = m_Phi*ts_in.m_A - m_R; % matrix containing gradients

        %for n = 1:N, for n_prime = 1:N
                % af_nnp[t] = a_nnp[t] - alpha_t*v_nnp[t]
                % compute a_nnp[t+1] via (24)
        % end
               
        % update buffer (shift buffer and add v_y at the end)
        ts_out = ts_in.updateBuffer(v_y);
        ts_out.m_Phi = m_Phi;
        ts_out.m_R   = m_R;
                
        alpha = obj.h_stepsize(ts_out);
        m_Af = ts_in.m_A - alpha*m_V;
        [ts_out.m_A, v_norms] = ...
            obj.group_soft_thresholding(m_Af, alpha);
        
        if nargout > 2
            m_B = obj.calculate_m_B(m_Af, v_norms, alpha);
        end
    end
    
    function [m_A_out, v_norms] = group_soft_thresholding(obj, m_Af, alpha)
        % group_soft_thresholding performs group- soft thresholding (just
        % like in Group Lasso) over each row of the m_A matrix
        % this is a prox operator
        % Input:  m_Af PN x N matrix forward iterate of Tirso
        %         alpha: step size that was used in the forward iteration
        % Output: m_A_out: PN x N matrix where each row is the result of 
        %                  the group- soft thresholding
        %         v_norms: N*N -vector containing the norms of each FIR
        %                  impulse response
        assert(isscalar(alpha), 'not ready for nonscalar step size')
        P = obj.order;
        N = obj.noOfNodes;
        
        m_lambdas = obj.get_m_lambdas;
        
        %% original implementation
        m_Af_fat = reshape(m_Af, [P N*N]);
        v_norms = sqrt(sum(m_Af_fat.^2));
        v_factors = max(0, 1-(alpha*m_lambdas(:)')./v_norms);
        % multiply each column of Af by its corresponding factor:
        m_Aout_fat = m_Af_fat.*repmat(v_factors, [P, 1]); 
        m_A_out = reshape(m_Aout_fat, [P*N N]);

        %% alternative implementation 
        t_Af = reshape(m_Af, [P N N]);
        m_norms = squeeze(vecnorm(t_Af));
        m_factors = max(0, 1-(alpha*m_lambdas)./m_norms);
%         assert (norm(m_factors(:)-v_factors(:))<1e-8)
        t_Aout = t_Af.*reshape(m_factors, [1, N, N]);
        m_A_out2 = reshape(t_Aout, [P*N N]);
        
%         assert(norm(vec(m_A_out2 - m_A_out)) < 1e-8);
        
    end
    
    function m_lambdas_out = get_m_lambdas(obj)
        % get_m_lambdas produces a matrix containing the regularization
        % Output: m_lambdas N x N matrix
        N = obj.noOfNodes;
                
        if isscalar(obj.regPar)
            m_lambdas = obj.regPar*(ones(N));
        elseif iscolumn(obj.regPar)
            m_lambdas = repmat(obj.regPar, [1 N]);
        elseif isrow(obj.regPar)
            warning('regPar should be a column vector, not a row. I transpose it for you.')
            m_lambdas = repmat(obj.regPar', [1 N]);
        elseif ismatrix(obj.regPar)        
            error('Matrix reg parameters not implemented yet')
            % TODO: more advanced version (one different lambda per link)
        else, error('invalid format of regPar property')
        end

        
        if obj.b_shrinkSelfLoops % Group Lasso that does not differentiate 
            % between inter-node links and self loops
            m_lambdas_out = m_lambdas;
        else % Bolstad, i.e., self loops do not suffer shrinking
            m_lambdas_out = m_lambdas.*(1-eye(N));
        end

    end
end
    
end