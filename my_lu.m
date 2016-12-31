function [p, lu_mat, q] = my_lu(A, is_full)
	n = rows(A);
	if n <= 1
		lu_mat = A;
		p = ones(1,1);
		q = ones(1,1);
	else
		#find max
		[max_rows, max_rows_arg] = max(A(1:n,1:1));
		[max_cols, max_cols_arg] = max(A(1:1, 1:n));
		p = eye(n);
		q = eye(n);
		if max_rows > max_cols || !is_full
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
		lu_mat(1,1:n) = A(1, 1:n);
		l = A(2:n,1) / A(1,1);
		lu_mat(2:n, 1) = l;
		A_sub = A(2:n, 2:n) - l * A(1,2:n);
		[sub_p, lu_mat(2:n,2:n), sub_q] = my_lu(A_sub, is_full);


		# transpose calculated rows and columns
		lu_mat(2:n,1) = sub_p * lu_mat(2:n, 1);
		lu_mat(1,2:n) = lu_mat(1, 2:n) * sub_q;
		new_p = eye(n);
		new_p(2:n,2:n) = sub_p;
		p = p * new_p;
		new_q = eye(n);
		new_q(2:n,2:n) = sub_q;
		q = new_q * q;
	endif
endfunction

function x = solve_low(lu_mat, b)
	n = rows(b);
	x = [b(1)];
	for i = 2:n
		x = [x; b(i) - lu_mat(i, 1:i-1) * x];
	endfor
endfunction


function x = solve_up(lu_mat, b)
	n = rows(b);
	x = [b(n) / lu_mat(n,n)];
	for i = n - 1:-1:1
		x = [(b(i) - lu_mat(i, i + 1:n) * x)/lu_mat(i,i) ;x];
	endfor
endfunction


function x = lu_solve(lu_mat, b)
	#solve low left
	z = solve_low(lu_mat, b);
	x = solve_up(lu_mat, z);
endfunction


function x = solve(A, b)
	[p, lu_mat, q] = my_lu (A, false);
	y = lu_solve(lu_mat, p * b);
	x = q * y;
endfunction


function [A, p, l, u, q ] = solve_6()
	n = 6;
	A = zeros(n,n);
	for i = 1:n
		for j = 1:n
			A(i,j) = 1 / (i + j -1);
		endfor
	endfor
	[p,lu_mat,q]=my_lu(A, false);
	l = tril(lu_mat, -1) + eye(n);
	u = triu(lu_mat,0);
endfunction


function easy_check()
	#generate data
	display("Generating data ...");
	A = [[1 4]; [2 3]];
	b = [11; 12];
	display("A = ");
	display(A);
	display("b = ");
	display(b);

	display("Solving ...");
	x = solve(A, b);

	display("x = ");
	display(x);
endfunction


function check_6(is_full)
	display("Check with A_6 = ");
	if is_full
		display("Full pivoting ...");
	else
		display("Partial pivoting ...");
	endif
	[A, p, l, u,q ] = solve_6();
	display("p * A * q = l * u");
	display(p);
	display("*");
	display(A);
	display("*");
	display(q);
	display(" = ");
	display(l);
	display("*");
	display(u);
	display("<=>");
	display(p * A * q);
	display("=");
	display(l * u);
endfunction


function main()
	easy_check();
	check_6(false);
	check_6(true);
	display("Done!");
endfunction

main()
