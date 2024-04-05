//
// Created by jordan on 06/06/23.
//

#include "nautilus-refine.h"


Optimiser_simplex::Optimiser_simplex( double tolerance, int max_cycles, TYPE type, bool debug )
{
    tolerance_ = tolerance;
    max_cycles_ = max_cycles;
    type_ = type;
    debug_mode = debug;
}

std::vector<double> Optimiser_simplex::operator() ( const Target_fn_order_zero&
target_fn, const std::vector<std::vector<double> >& args ) const
{
    if (debug_mode)
        std::cout << "() operator called " << std::endl;
    enum STEP { UNKN, EXTN, NRML, CTRN, CTRX };

    // check input
    int size = target_fn.num_params();
    if ( args.size() != size+1 )
        clipper::Message::message( clipper::Message_fatal( "Optimiser_simplex: parameter size mismatch" ) );
    for ( int i = 0; i < args.size(); i++ )
        if ( args[i].size() != size )
            clipper::Message::message( clipper::Message_fatal( "Optimiser_simplex: parameter size mismatch" ) );

    // make arrays
    std::vector<std::vector<double> > params( size + 1 );
    std::vector<double> fn( size + 1 );
    // calc target fn
    for ( int i = 0; i < args.size(); i++ ) {
        params[i] = args[i];
        fn[i] = target_fn( args[i] );
    }

    // simplex loop
    int worst_index, best_index;
    worst_index = best_index = 0;
    double f0(0.0), f1(0.0), f2(0.0);
    for ( int cyc = 0; cyc < max_cycles_; cyc++ ) {
        if (debug_mode)
            std::cout << cyc << "/" << max_cycles_ << std::endl;

        if ( debug_mode )  // DEBUG OUTPUT
            for ( int i = 0; i < params.size(); i++ ) {
                std::cout << "Index " << i << "  score of i " << clipper::String( fn[i] ) << "  ";
                for ( int j = 0; j < size; j++ )
                    std::cout << "params[" << i << "][" << j << "] " << clipper::String( params[i][j] ) << "\t";
                std::cout << "\n";
            }
        // find worst point
        worst_index = best_index = 0;
        for ( int i = 0; i < params.size(); i++ ) {
            if ( fn[i] > fn[worst_index] ) worst_index = i;
            if ( fn[i] < fn[best_index] ) best_index = i;
        }

        if (debug_mode)
            std::cout << "worst " << worst_index << "\tbest " << best_index << std::endl;

        // termination condition
        if (fn[worst_index] - fn[best_index] < tolerance_ )
        {
            if (debug_mode)
                std::cout << "Termination hit " << std::endl;
            break;
        }
        // find centroid of the rest
        std::vector<double> centr( size, 0.0 ),
                shift( size ), t0( size ), t1( size ), t2( size );
        double sumweight = 0.0;

        for ( int i = 0; i < params.size(); i++ )
            if (i != worst_index ) {
                double weight = 1.0;
                if ( type_ == GRADIENT ) {
                    double r2 = 0.0;
                    for ( int j = 0; j < size; j++ )
                        r2 += clipper::Util::sqr( params[i][j] - params[worst_index][j] );
                    weight = (fn[worst_index] - fn[i] ) / sqrt(r2 );
                }
                for ( int j = 0; j < size; j++ )
                    centr[j] += weight * params[i][j];
                sumweight += weight;
            }

        for ( int j = 0; j < size; j++ ) centr[j] = centr[j] / sumweight;
        for ( int j = 0; j < size; j++ ) shift[j] = centr[j] - params[worst_index][j];
        // calculate first trial point
        for ( int j = 0; j < size; j++ ) {
//            std::cout << centr[j] << "\t" << shift[j] << std::endl;
            t0[j] = centr[j] - 0.5 * shift[j];
            t1[j] = centr[j] + 1.0 * shift[j];
            t2[j] = centr[j] + 2.0 * shift[j];
        }


//        std::cout << "Centr " << centr[0] << " " << centr[1] << " " << centr[2] << std::endl;
//        std::cout << "Shift " << shift[0] << " " << shift[1] << " " << shift[2] << std::endl;

        f1 = target_fn( t1 );
        // simplex conditions
        STEP step = UNKN;

        if ( !clipper::Util::is_nan( f1 ) ) { // new point is valid
            if ( f1 < fn[worst_index] ) {    // new point is better than worst
                if ( f1 < fn[best_index] ) {  // new point is better than best
                    f2 = target_fn( t2 );
                    if ( !clipper::Util::is_nan( f2 ) ) {
                        if ( f2 < f1 ) {  // extended step is best
                            step = EXTN;
                        } else {  // normal step better than extended step
                            step = NRML;
                        }
                    } else {  // extended step is invalid
                        step = NRML;
                    }
                } else {  // normal step is neither best nor worst (default)
                    step = NRML;
                }
            }  // else normal step is worse
        }  // else normal step is invalid
        if ( step == UNKN ) {    // normal step is worse or invalid
            f0 = target_fn( t0 );  // contraction step (always valid)
            if ( f0 < fn[worst_index] ) {   // contraction step is better than worst
                step = CTRN;
            } else {               // contraction step is worse
                step = CTRX;
            }
        }
        if      ( step == EXTN ) { params[worst_index] = t2; fn[worst_index] = f2; }  // extension
        else if ( step == NRML ) { params[worst_index] = t1; fn[worst_index] = f1; }  // normal
        else if ( step == CTRN ) { params[worst_index] = t0; fn[worst_index] = f0; }  // contraction
        else {  // otherwise contract all towards best (slow)
            for ( int i = 0; i < params.size(); i++ )
                if (i != best_index ) {
                    for ( int j = 0; j < size; j++ )
                        params[i][j] = 0.5 * ( params[i][j] + params[best_index][j] );
                    fn[i] = target_fn( params[i] );
                }
        }
        if ( debug_mode ) {  // DEBUG OUTPUT
            if      ( step == EXTN ) std::cout << "Extn-step\n";
            else if ( step == NRML ) std::cout << "Nrml-step\n";
            else if ( step == CTRN ) std::cout << "Ctrn-step\n";
            else                     std::cout << "Ctrx-step\n";
        }
    }
    return params[best_index];
}


double Target_fn_refine_phosphate::operator()(const std::vector<double> &args) const {
    clipper::Coord_orth coord = {args[0], args[1], args[2]};
    return -m_xmap->interp<clipper::Interp_cubic>(coord.coord_frac(m_xmap->cell()));

}

clipper::Coord_orth Target_fn_refine_phosphate::refine(const clipper::Coord_orth &coord) {

    std::vector<double> args = {coord.x(), coord.y(), coord.z()};

    std::vector<std::vector<double>> args_init;

    args_init.push_back( args );
    args_init.push_back({args[0]+m_step, args[1], args[2]});
    args_init.push_back({args[0], args[1]+m_step, args[2]});
    args_init.push_back({args[0], args[1], args[2]+m_step});

    double tol = 0.005 * (*this)( args_init[0] );
    Optimiser_simplex optimiser_simplex(tol, 100, Optimiser_simplex::NORMAL, false);
    std::vector<double> args_refined = optimiser_simplex(*this, args_init);
    return {args_refined[0], args_refined[1], args_refined[2]};

}

double Target_fn_refine_fragment::operator()(const std::vector<double> &args) const {

    clipper::Euler_ccp4 euler = {args[0], args[1], args[2]};
    clipper::Rotation rotation(euler);
    clipper::Mat33<> rot_mat = rotation.matrix();

    clipper::RTop_orth rtop = {rot_mat, m_translation};

    NucleicAcidDB::ChainFull test_fragment = m_chain;
    test_fragment.transform(rtop);

    float score = score_fragment(test_fragment, *m_xmap);
    return -score;
}


//double Target_fn_refine_fragment::operator()(const std::vector<double> &args) const {
//
//    clipper::Euler_ccp4 euler = {args[0], args[1], args[2]};
//    clipper::Rotation rotation(euler);
//    clipper::Mat33<> rot_mat = rotation.matrix();
//
//    clipper::RTop_orth rtop = {rot_mat, m_translation};
//
//    NucleicAcidDB::ChainFull test_fragment = m_chain;
//    test_fragment.transform(rtop);
//
//    int total_count = 0;
//    float total_score = 0.0f;
//    for (int i = 0; i < test_fragment.size(); i++) {
//
//        float score = 0.0f;
//        auto chain = test_fragment[i];
//
//        if (!chain.O5p1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.O5p1.coord_frac(m_xmap->cell()));
//        if (!chain.C5p1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.C5p1.coord_frac(m_xmap->cell()));
//        if (!chain.C4p1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.C4p1.coord_frac(m_xmap->cell()));
//        if (!chain.O4p1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.O4p1.coord_frac(m_xmap->cell()));
//        if (!chain.C3p1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.C3p1.coord_frac(m_xmap->cell()));
//        if (!chain.O3p1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.O3p1.coord_frac(m_xmap->cell()));
//        if (!chain.C2p1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.C2p1.coord_frac(m_xmap->cell()));
//        if (!chain.C1p1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.C1p1.coord_frac(m_xmap->cell()));
//        if (!chain.C2_1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.C2_1.coord_frac(m_xmap->cell()));
//        if (!chain.C4_1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.C4_1.coord_frac(m_xmap->cell()));
//        if (!chain.C5_1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.C5_1.coord_frac(m_xmap->cell()));
//        if (!chain.C6_1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.C6_1.coord_frac(m_xmap->cell()));
//        if (!chain.C8_1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.C8_1.coord_frac(m_xmap->cell()));
//        if (!chain.N1_1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.N1_1.coord_frac(m_xmap->cell()));
//        if (!chain.N2_1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.N2_1.coord_frac(m_xmap->cell()));
//        if (!chain.N3_1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.N3_1.coord_frac(m_xmap->cell()));
//        if (!chain.N4_1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.N4_1.coord_frac(m_xmap->cell()));
//        if (!chain.N6_1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.N6_1.coord_frac(m_xmap->cell()));
//        if (!chain.N7_1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.N7_1.coord_frac(m_xmap->cell()));
//        if (!chain.N9_1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.N9_1.coord_frac(m_xmap->cell()));
//        if (!chain.O2_1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.O2_1.coord_frac(m_xmap->cell()));
//        if (!chain.O6_1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.O6_1.coord_frac(m_xmap->cell()));
//
//        test_fragment[i].score = score;
//        total_score += score;
//        total_count += chain.get_mmonomer().size();
//    }
////    std::cout << euler.format() << "-> " << -total_score/total_count << std::endl;
//
//    return -total_score/total_count;
//}

clipper::RTop_orth
Target_fn_refine_fragment::refine() {

//    clipper::Mat33<> mat = rtop.rot();
//    clipper::Rotation rotation(mat);
//    clipper::Euler_ccp4 euler = rotation.euler_ccp4();
//    std::vector<double> args = {euler.alpha(), euler.beta(), euler.gamma()};
    std::vector<double> args = {0, 0, 0};

    std::vector<std::vector<double>> args_init;

    float step = clipper::Util::d2rad(2);

    args_init.push_back( args );
    args_init.push_back({args[0]+step, args[1], args[2]});
    args_init.push_back({args[0], args[1]+step, args[2]});
    args_init.push_back({args[0], args[1], args[2]+step});

    double tol = 0.0005 * (*this)( args_init[0] );
    Optimiser_simplex optimiser_simplex(tol, 20, Optimiser_simplex::NORMAL, false );
    std::vector<double> args_refined = optimiser_simplex(*this, args_init);

    clipper::Euler_ccp4 refined_euler = {args_refined[0], args_refined[1], args_refined[2]};
    clipper::Rotation refined_rotation(refined_euler);
    clipper::Mat33<> refined_mat = refined_rotation.matrix();

    return {refined_mat, m_translation};
}



double Target_fn_refine_fragment_trn::operator()(const std::vector<double> &args) const {

    clipper::Vec3<> new_trn = {args[0], args[1], args[2]};

    clipper::RTop_orth rtop = {clipper::Mat33<>::identity(), new_trn};

    NucleicAcidDB::ChainFull test_fragment = m_chain;
    test_fragment.transform(rtop);

    int total_count = 0;
    float total_score = 0.0f;
    for (int i = 0; i < test_fragment.size(); i++) {

        float score = 0.0f;
        auto chain = test_fragment[i];

        if (!chain.O5p1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.O5p1.coord_frac(m_xmap->cell()));
        if (!chain.C5p1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.C5p1.coord_frac(m_xmap->cell()));
        if (!chain.C4p1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.C4p1.coord_frac(m_xmap->cell()));
        if (!chain.O4p1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.O4p1.coord_frac(m_xmap->cell()));
        if (!chain.C3p1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.C3p1.coord_frac(m_xmap->cell()));
        if (!chain.O3p1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.O3p1.coord_frac(m_xmap->cell()));
        if (!chain.C2p1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.C2p1.coord_frac(m_xmap->cell()));
        if (!chain.C1p1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.C1p1.coord_frac(m_xmap->cell()));
        if (!chain.C2_1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.C2_1.coord_frac(m_xmap->cell()));
        if (!chain.C4_1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.C4_1.coord_frac(m_xmap->cell()));
        if (!chain.C5_1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.C5_1.coord_frac(m_xmap->cell()));
        if (!chain.C6_1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.C6_1.coord_frac(m_xmap->cell()));
        if (!chain.C8_1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.C8_1.coord_frac(m_xmap->cell()));
        if (!chain.N1_1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.N1_1.coord_frac(m_xmap->cell()));
        if (!chain.N2_1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.N2_1.coord_frac(m_xmap->cell()));
        if (!chain.N3_1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.N3_1.coord_frac(m_xmap->cell()));
        if (!chain.N4_1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.N4_1.coord_frac(m_xmap->cell()));
        if (!chain.N6_1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.N6_1.coord_frac(m_xmap->cell()));
        if (!chain.N7_1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.N7_1.coord_frac(m_xmap->cell()));
        if (!chain.N9_1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.N9_1.coord_frac(m_xmap->cell()));
        if (!chain.O2_1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.O2_1.coord_frac(m_xmap->cell()));
        if (!chain.O6_1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.O6_1.coord_frac(m_xmap->cell()));

        test_fragment[i].score = score;
        total_score += score;
        total_count += chain.get_mmonomer().size();
    }

    std::cout << new_trn.format() << "\t" << -total_score/total_count << "\t" << test_fragment.size() << std::endl;

    return -total_score/total_count;
}

clipper::RTop_orth
Target_fn_refine_fragment_trn::refine() {

    std::vector<double> args = {m_translation[0], m_translation[1], m_translation[2]};

    std::vector<std::vector<double>> args_init;

    float step = 0.2;

    args_init.push_back( args );
    args_init.push_back({args[0]+step, args[1], args[2]});
    args_init.push_back({args[0], args[1]+step, args[2]});
    args_init.push_back({args[0], args[1], args[2]+step});

    double tol = 0.0005 * (*this)( args_init[0] );
    Optimiser_simplex optimiser_simplex(tol, 100, Optimiser_simplex::GRADIENT, true );
    std::vector<double> args_refined = optimiser_simplex(*this, args_init);

    clipper::Vec3<> refined_trn = {args_refined[0], args_refined[1], args_refined[2]};

//    std::cout << m_translation.format() << "\t" << refined_trn.format() << std::endl;

    return {clipper::Mat33<>::identity(), refined_trn};

}

// REFINE FRAGMENT COORDINATES

double RefineFragmentCoordinates::operator()(const std::vector<double> &args) const {

    clipper::Vec3<> new_translation = {args[3], args[4], args[5]};

    clipper::Euler_ccp4 refined_euler = {args[0], args[1], args[2]};
    clipper::Rotation refined_rotation(refined_euler);
    clipper::Mat33<> refined_rotation_matrix = refined_rotation.matrix();

    clipper::RTop_orth rtop = {refined_rotation_matrix, new_translation};

    NucleicAcidDB::ChainFull test_fragment = m_chain;
    test_fragment.transform(rtop);

    int total_count = 0;
    float total_score = 0.0f;
    for (int i = 0; i < test_fragment.size(); i++) {

        float score = 0.0f;
        auto chain = test_fragment[i];

        if (!chain.O5p1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.O5p1.coord_frac(m_xmap->cell()));
        if (!chain.C5p1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.C5p1.coord_frac(m_xmap->cell()));
        if (!chain.C4p1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.C4p1.coord_frac(m_xmap->cell()));
        if (!chain.O4p1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.O4p1.coord_frac(m_xmap->cell()));
        if (!chain.C3p1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.C3p1.coord_frac(m_xmap->cell()));
        if (!chain.O3p1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.O3p1.coord_frac(m_xmap->cell()));
        if (!chain.C2p1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.C2p1.coord_frac(m_xmap->cell()));
        if (!chain.C1p1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.C1p1.coord_frac(m_xmap->cell()));
        if (!chain.C2_1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.C2_1.coord_frac(m_xmap->cell()));
        if (!chain.C4_1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.C4_1.coord_frac(m_xmap->cell()));
        if (!chain.C5_1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.C5_1.coord_frac(m_xmap->cell()));
        if (!chain.C6_1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.C6_1.coord_frac(m_xmap->cell()));
        if (!chain.C8_1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.C8_1.coord_frac(m_xmap->cell()));
        if (!chain.N1_1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.N1_1.coord_frac(m_xmap->cell()));
        if (!chain.N2_1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.N2_1.coord_frac(m_xmap->cell()));
        if (!chain.N3_1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.N3_1.coord_frac(m_xmap->cell()));
        if (!chain.N4_1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.N4_1.coord_frac(m_xmap->cell()));
        if (!chain.N6_1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.N6_1.coord_frac(m_xmap->cell()));
        if (!chain.N7_1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.N7_1.coord_frac(m_xmap->cell()));
        if (!chain.N9_1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.N9_1.coord_frac(m_xmap->cell()));
        if (!chain.O2_1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.O2_1.coord_frac(m_xmap->cell()));
        if (!chain.O6_1.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.O6_1.coord_frac(m_xmap->cell()));
        if (!chain.P.is_null()) score += m_xmap->interp<clipper::Interp_cubic>(chain.P.coord_frac(m_xmap->cell()));
        if (!chain.P.is_null()) {
            float density_predicted_score = m_phosphate_map->interp<clipper::Interp_cubic>(chain.P.coord_frac(m_xmap->cell()));
            if (density_predicted_score < 0.1) {
                score -= 100;
            }
        }

        test_fragment[i].score = score;
        total_score += score;
        total_count += chain.get_mmonomer().size();
    }

    //    std::cout << new_translation.format() << "\t" << refined_euler.format()<< std::endl;

    return -total_score/total_count;
}

clipper::RTop_orth RefineFragmentCoordinates::refine() {
    std::vector<double> args = {0, 0, 0, m_translation[0], m_translation[1], m_translation[2]};

    std::vector<std::vector<double>> args_init;

    args_init.push_back( args );
    args_init.push_back({args[0]+m_rotation_step, args[1], args[2], args[3], args[4], args[5] });
    args_init.push_back({args[0], args[1]+m_rotation_step, args[2], args[3], args[4], args[5] });
    args_init.push_back({args[0], args[1], args[2]+m_rotation_step, args[3], args[4], args[5] });
    args_init.push_back({args[0], args[1], args[2], args[3]+m_translation_step, args[4], args[5] });
    args_init.push_back({args[0], args[1], args[2], args[3], args[4]+m_translation_step, args[5] });
    args_init.push_back({args[0], args[1], args[2], args[3], args[4], args[5]+m_translation_step });

    double tol = 0.0005 * (*this)( args_init[0] );
    Optimiser_simplex optimiser_simplex(tol, 20, Optimiser_simplex::NORMAL, false );
    std::vector<double> args_refined = optimiser_simplex(*this, args_init);

    clipper::Euler_ccp4 refined_euler = {args_refined[0], args_refined[1], args_refined[2]};
    clipper::Rotation refined_rotation(refined_euler);
    clipper::Mat33<> refined_rotation_matrix = refined_rotation.matrix();

    clipper::Vec3<> refined_translation = {args_refined[3], args_refined[4], args_refined[5]};

    return {refined_rotation_matrix, refined_translation};

}

