function opt = generateRotFrameFrequencyVector( opt )

if isfield( opt, 'dt' )
    dt = opt.dt; % get spacing of points for optimization
else
    dt = 5e-6;
end
dt_tol = dt / 10;

if isfield( opt, 'tvec' )
    tvec = opt.tvec;
elseif isfield( opt, 'init_tvec' )
    tvec = opt.init_tvec;
end

if isfield( opt, 'constantRotatingFrame' )
    if opt.constantRotatingFrame
        opt.dwxyvec = opt.dwxy * ones( 1, opt.numTimePoints );
    else
        if ~isfield( opt, 'dwxyvec' )
            opt.dwxyvec = zeros( 1, opt.numTimePoints );

            % determine if there are any gaps in the timing array
            numPeriods = ( size( opt.timing_dwxy, 1 ) );
            pp = 1;
            while pp <= numPeriods

                if pp == 1
                    if opt.timing_dwxy( pp, 1 ) ~= 0
                        opt.timing_dwxy( pp, 1 ) = 0;
                    end
                end
                if pp == size( opt.timing_dwxy, 1 )
                    if opt.timing_dwxy( pp, 2 ) ~= opt.pulseLength
                        opt.timing_dwxy( pp, 2 ) = opt.pulseLength;
                    end
                end

                if pp < ( size( opt.timing_dwxy, 1 ) )
                    if abs( opt.timing_dwxy( pp, 2 ) - opt.timing_dwxy( pp+1, 1 ) ) > dt_tol

                        opt.timing_dwxy = [...
                            opt.timing_dwxy( 1:pp, : );...
                            [ opt.timing_dwxy( pp, 2 ), opt.timing_dwxy( pp+1, 1 ) ];...
                            opt.timing_dwxy( (pp+1):end, : );...
                            ];

                        opt.vals_dwxy = [...
                            opt.vals_dwxy( 1:pp, : );...
                            opt.vals_dwxy( pp );...
                            opt.vals_dwxy( (pp+1):end, : );...
                            ];
                    end
                end

                pp = pp + 1;
                numPeriods = ( size( opt.timing_dwxy, 1 ) );

            end

            for pp = 1:size( opt.timing_dwxy, 1 )
                perIdx = ( ( tvec ) > ( opt.timing_dwxy( pp, 1 ) - dt_tol ) ) &...
                    ( ( tvec ) < ( opt.timing_dwxy( pp, 2 ) + dt_tol ) );
                opt.dwxyvec( perIdx ) = opt.vals_dwxy( pp );
            end
        end

    end
else
    opt.dwxyvec = zeros( 1, opt.numTimePoints );
    opt.constantRotatingFrame = true;
end

end