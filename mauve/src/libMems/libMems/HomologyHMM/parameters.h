#ifndef __HomologyHMM_parameters_h__
#define __HomologyHMM_parameters_h__

#include "homology.h"

Params getHoxdParams();
Params getAdaptedHoxdMatrixParameters( double gc_content );
void adaptToPercentIdentity( Params& params, double pct_identity );

inline
Params getHoxdParams()
{
    static Params* params = NULL;
	if( params == NULL )
	{
		params = new Params();
		params->iStartHomologous = 0.5;
		params->iGoHomologous = 0.00001;
		params->iGoUnrelated = 0.0000001;
		params->iGoStopFromUnrelated = 0.00000001;
		params->iGoStopFromHomologous = 0.00000001;

		// original values from Chiaromonte et al supplied by Webb Miler
		params->aEmitHomologous[0] = 0.1723*2;		//a:a, t:t
		params->aEmitHomologous[1] = 0.1462*2;	//c:c, g:g
		params->aEmitHomologous[2] = 0.0180*4;	//a:c, c:a, g:t, t:g
		params->aEmitHomologous[3] = 0.0426*4;	//a:g, g:a, c:t, t:c
		params->aEmitHomologous[4] = 0.0186*2;	//a:t, t:a
		params->aEmitHomologous[5] = 0.0142*2;	//g:c, c:g
		params->aEmitHomologous[6] = 0.004461;	// gap open (from an e. coli y pestis alignment)
		// gap extend // 0.050733
		params->aEmitHomologous[7] = 1 - (params->aEmitHomologous[0] + params->aEmitHomologous[1] + params->aEmitHomologous[2] +
				params->aEmitHomologous[3] + params->aEmitHomologous[4] + params->aEmitHomologous[5] + params->aEmitHomologous[6]);


		params->aEmitUnrelated[0] = 0.12818742714404662781015820149872;	// a:a, t:t
		params->aEmitUnrelated[1] = 0.10493347210657785179017485428807;	// c:c, g:g
		params->aEmitUnrelated[2] = 0.11597910074937552039966694421313;	// a:c, c:a
		params->aEmitUnrelated[3] = params->aEmitUnrelated[2];
		params->aEmitUnrelated[4] = params->aEmitUnrelated[0];
		params->aEmitUnrelated[5] = params->aEmitUnrelated[1]; 
		params->aEmitUnrelated[6] = 0.0483;	// gap open (derived by aligning a 48%GC sequence with 
										// its reverse--not complement--to derive expected gap frequencies in 
										// unrelated sequence)
		// gap extend // 0.2535
		params->aEmitUnrelated[7] = 1 - (params->aEmitUnrelated[0] + params->aEmitUnrelated[1] + params->aEmitUnrelated[2] +
				params->aEmitUnrelated[3] + params->aEmitUnrelated[4] + params->aEmitUnrelated[5] + params->aEmitUnrelated[6]);
	}

	return *params;
}


/**
 * Adapts an emission matrix to an arbitrary nucleotide composition
 * @param gc_content	The fraction of the genome which is G/C
 */
inline
Params getAdaptedHoxdMatrixParameters( double gc_content )
{
	Params params;
    double s = 0.03028173853;
    double at_content = 1-gc_content;
    double norm_factor = 0.0;

	double gO_unrelated = 0.0483;
	double gE_unrelated = 0.2535;

	double gO_homologous = 0.004461;
	double gE_homologous = 0.050733;

    // Unrelated state emission probabilities
    // use AT/GC background frequency instead of mononucleotide frequency since that is how it is described in the manuscript
    params.aEmitUnrelated[0] = (at_content/2)*(at_content/2)+(at_content/2)*(at_content/2); // a:a, t:t
    params.aEmitUnrelated[1] = (gc_content/2)*(gc_content/2)+(gc_content/2)*(gc_content/2); // c:c, g:g
    params.aEmitUnrelated[2] = (at_content/2)*(gc_content/2)+(gc_content/2)*(at_content/2); //a:c, c:a, g:t, t:g
    params.aEmitUnrelated[3] = params.aEmitUnrelated[2]; //a:g, g:a, c:t, t:c
    params.aEmitUnrelated[4] = params.aEmitUnrelated[0]; //a:t, t:a 
    params.aEmitUnrelated[5] = params.aEmitUnrelated[1]; //g:c, c:g 
    
    
    norm_factor = (1-(gO_unrelated+gE_unrelated))/(params.aEmitUnrelated[0] + params.aEmitUnrelated[1] +params.aEmitUnrelated[2] + params.aEmitUnrelated[3] 
                        + params.aEmitUnrelated[4] + params.aEmitUnrelated[5] );

    //NORMALIZE the values
    params.aEmitUnrelated[0] = params.aEmitUnrelated[0]*norm_factor;
    params.aEmitUnrelated[1] = params.aEmitUnrelated[1]*norm_factor;
    params.aEmitUnrelated[2] = params.aEmitUnrelated[2]*norm_factor;
    params.aEmitUnrelated[3] = params.aEmitUnrelated[3]*norm_factor;
    params.aEmitUnrelated[4] = params.aEmitUnrelated[4]*norm_factor;
    params.aEmitUnrelated[5] = params.aEmitUnrelated[5]*norm_factor;
    params.aEmitUnrelated[6] = gO_unrelated;// gap open 
    params.aEmitUnrelated[7] = 1 - (params.aEmitUnrelated[0] + params.aEmitUnrelated[1] + params.aEmitUnrelated[2] + params.aEmitUnrelated[3] 
                        + params.aEmitUnrelated[4] + params.aEmitUnrelated[5] + params.aEmitUnrelated[6]);

    //USE PRE-NORMALIZED VALUES (from the HOXD matrix)!!
    double H_AA = 0.1723*2;		//a:a, t:t
    double H_CC = 0.1462*2;		//c:c, g:g
    double H_AC = 0.0180*4;	//a:c, c:a, g:t, t:g
    double H_AG = 0.0426*4;	//a:g, g:a, c:t, t:c
    double H_AT = 0.0186*2;	//a:t, t:a
    double H_CG = 0.0142*2;	//g:c, c:g

    // Homologous state emission probabilities 
    params.aEmitHomologous[0] = (at_content/0.525)*H_AA; // a:a, t:t
    params.aEmitHomologous[1] = (gc_content/0.475)*H_CC; // c:c, g:g
    params.aEmitHomologous[2] = H_AC; //a:c, c:a, g:t, t:g
    params.aEmitHomologous[3] = H_AG; //a:g, g:a, c:t, t:c
    params.aEmitHomologous[4] = (at_content/0.525)*H_AT; //a:t, t:a 
    params.aEmitHomologous[5] = (gc_content/0.475)*H_CG; //g:c, c:g 

    
    norm_factor = (1-(gO_homologous+gE_homologous))/(params.aEmitHomologous[0] + params.aEmitHomologous[1] + params.aEmitHomologous[2] + params.aEmitHomologous[3] 
                    + params.aEmitHomologous[4] + params.aEmitHomologous[5]);
    
    //NORMALIZE the values
    params.aEmitHomologous[0] = params.aEmitHomologous[0]*norm_factor;
    params.aEmitHomologous[1] = params.aEmitHomologous[1]*norm_factor;
    params.aEmitHomologous[2] = params.aEmitHomologous[2]*norm_factor;
    params.aEmitHomologous[3] = params.aEmitHomologous[3]*norm_factor;
    params.aEmitHomologous[4] = params.aEmitHomologous[4]*norm_factor;
    params.aEmitHomologous[5] = params.aEmitHomologous[5]*norm_factor;
    params.aEmitHomologous[6] = gO_homologous;// gap open
    params.aEmitHomologous[7] = 1 - (params.aEmitHomologous[0] + params.aEmitHomologous[1] + params.aEmitHomologous[2] + params.aEmitHomologous[3] 
                        + params.aEmitHomologous[4] + params.aEmitHomologous[5] + params.aEmitHomologous[6]);


	// set state transition probabilities
	params.iStartHomologous = 0.5;
	params.iGoHomologous = 0.00001;
	params.iGoUnrelated = 0.0000001;

	params.iGoStopFromHomologous = 0.0000001;
	params.iGoStopFromUnrelated = 0.0000001;

	return params;
}

inline
void adaptToPercentIdentity( Params& params, double pct_identity )
{
	if( pct_identity <= 0 || pct_identity > 1 )
		throw "Bad pct identity";		// error condition
	// normalize new pct identity to gap content
	double gapnorm_pct_id = pct_identity * (1.0 - params.aEmitHomologous[6] - params.aEmitHomologous[7]);
	// calculate the previous expected identity as H_AA + H_CC
	double prev_pct_id = params.aEmitHomologous[0] + params.aEmitHomologous[1];
	double id_diff = prev_pct_id - gapnorm_pct_id;
	// spread id_diff proportionally among other substitutions
	double rest_sum = params.aEmitHomologous[2] + params.aEmitHomologous[3] + 
		params.aEmitHomologous[4] + params.aEmitHomologous[5];
	params.aEmitHomologous[2] += id_diff * params.aEmitHomologous[2] / rest_sum;
	params.aEmitHomologous[3] += id_diff * params.aEmitHomologous[3] / rest_sum;
	params.aEmitHomologous[4] += id_diff * params.aEmitHomologous[4] / rest_sum;
	params.aEmitHomologous[5] += id_diff * params.aEmitHomologous[5] / rest_sum;

	params.aEmitHomologous[0] -= id_diff * params.aEmitHomologous[0] / prev_pct_id;
	params.aEmitHomologous[1] -= id_diff * params.aEmitHomologous[1] / prev_pct_id;
}

#endif	// __HomologyHMM_parameters_h__

