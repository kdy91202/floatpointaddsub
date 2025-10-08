module fp_addsub(in_a, in_b, sub, out);
    input [31:0] in_a, in_b;
    input sub;
    output [31:0] out;

    wire sign_a, sign_b;
    wire [7:0] exp_a, exp_b;
    wire [22:0] frac_a, frac_b;
    assign {sign_a, exp_a, frac_a} = in_a;
    assign {sign_b, exp_b, frac_b} = in_b;

    wire magnitude_a_gt_b = (exp_a > exp_b) || ((exp_a == exp_b) && (frac_a > frac_b));
    wire [23:0] sig_a = {1'b1, frac_a};
    wire [23:0] sig_b = {1'b1, frac_b};

    wire [7:0] exp_ge = magnitude_a_gt_b ? exp_a : exp_b;
    wire [7:0] exp_lt = magnitude_a_gt_b ? exp_b : exp_a;
    wire [23:0] sig_ge = magnitude_a_gt_b ? sig_a : sig_b;
    wire [23:0] sig_lt = magnitude_a_gt_b ? sig_b : sig_a;
    wire sign_ge = magnitude_a_gt_b ? sign_a : sign_b;
    wire [7:0] exp_diff = exp_ge - exp_lt;

    wire [23:0] sig_lt_shifted = sig_lt >> exp_diff;

    wire is_add = (~sub && ~(sign_a ^ sign_b)) | (sub & (sign_a ^ sign_b));
    wire [24:0] raw_mantissa = is_add ? ({1'b0, sig_ge} + {1'b0, sig_lt_shifted}) : ({1'b0, sig_ge} - {1'b0, sig_lt_shifted});
    wire raw_sign = sign_ge;
    wire [7:0] raw_exp = exp_ge;
    
    wire [7:0] exp_out;
    wire [22:0] frac_out;
    fp_norm u_norm (raw_exp, raw_mantissa, exp_out, frac_out);
    
    assign out = {raw_sign, exp_out, frac_out};
endmodule
	
module fp_norm(exp_in, mantissa_in, exp_out, mantissa_out);
    localparam M_WIDTH = 25;
    input signed [7:0] exp_in;
    input [M_WIDTH-1:0] mantissa_in;
    output signed [7:0] exp_out;
    output [22:0] mantissa_out;

    localparam EXT_WIDTH = 32;
    wire [EXT_WIDTH-1:0] mantissa_extended = {{(EXT_WIDTH - M_WIDTH){1'b0}}, mantissa_in};
    wire [$clog2(EXT_WIDTH)-1:0] lz_count_extended;
    wire valid;

    lod #(EXT_WIDTH) u_lod (mantissa_extended, lz_count_extended, valid);

    wire signed [$clog2(M_WIDTH):0] lz_count_actual = lz_count_extended - (EXT_WIDTH - M_WIDTH);
    
    wire signed [$clog2(M_WIDTH):0] shift_left = lz_count_actual - 1;
    assign exp_out = exp_in - shift_left;
    wire [M_WIDTH-1:0] mantissa_norm = mantissa_in << shift_left;
    
    assign mantissa_out = mantissa_norm[22:0];
endmodule

// leading one position detector
module lod(in, p, v);
	parameter N=32;
	localparam W = $clog2(N);
	input [N-1:0] in;
	output [W-1:0] p;		// position
	output v;				// valid
	
	generate
		if (N==2) begin
			assign v = in[1] | in[0];
			assign p[0] = ~in[1] & in[0];	// 01 -> 1, 00, 10, 11 -> 0
		end
		else begin
			wire [W-2:0] pL, pH;
			wire vL, vH;
			lod #(N>>1) u1 (in[(N>>1)-1:0], pL, vL);		// lower half input
			lod #(N>>1) u2 (in[N-1:(N>>1)], pH, vH);		// upper half input
			assign v = vH | vL;
			assign p = vH ? { 1'b0, pH} : {vL, pL};
		end
	endgenerate
endmodule