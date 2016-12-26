function [p, lu_mat, q] = lu(A)
	n = rows(A);
	if n <= 1
		lu_mat = A;
		p = ones(1,1);
		q = ones(1,1);
	else
		#find max
		[max_rows, max_rows_arg] = max(A(1:n,1:1));
		[max_cols, max_cols_arg] = max(A(1:1, 1:n));
		p = ones(n,n);
		q = ones(n,n);
		if max_rows > max_cols
			p(1,1) = 0;
			p(1,max_rows_arg) = 1;
			p(max_rows_arg, max_rows_arg) = 0;
			p(max_rows_arg, 1) = 1;
		else
			q(1,1) = 0;
			q(max_cols_arg, 1) = 1;
			q(max_cols_arg, max_cols_arg) = 0;
			q(1, max_cols_arg) = 1;
		endif
		A = p * A * q;


		lu_mat = zeros(n,n);
		lu_mat(1:1,1:n) = A(1:1, 1:n);
		l = A(2:n,1:1) / A(1,1);
		lu_mat(2:n, 1:1) = l;
		A_sub = A(2:n, 2:n) - l * A(1:1,2:n)';
		[sub_p, lu_mat(2:n,2:n), sub_q] = lu_mat(A_sub);


		# transpose calculated rows and columns
		lu_mat(2:n,1:1) = sub_p * lu_mat(2:n, 1:1);
		lu_mat(1:1,2:n) = lu_mat(1:1, 2:n) * sub_q;
		p(2:n,2:n) = sub_p;
		q(2:n,2:n) = sub_q;
	endif
endfunction

function x = solve_low(lu_mat, b)
	x = b - sum(tril(lu_mat, -1)')';
endfunction


function x = solve_up(lu_mat, b)
	x = (b - sum(triu(lu_mat, 1))) ./ diag(lu_mat);
endfunction


function x = lu_solve(lu_mat, b)
	#solve low left
	z = solve_low(lu_mat, b);
	x = solve_up(lu_mat, z);
endfunction


function x = solve(A, b)
	[p, lu_mat, q] = lu(A);
	y = lu_solve(lu_mat, P * b);
	x = q * y;
endfunction




function main()
	#generate data
	display("Generating data ...");
	A = [[1 4]; [2 3]];
	b = [11; 12];
	#display(A);




	display("Solving ...");
	x = solve(A, b);


	display("x = ", x);
	display("Done!");
endfunction




main()
