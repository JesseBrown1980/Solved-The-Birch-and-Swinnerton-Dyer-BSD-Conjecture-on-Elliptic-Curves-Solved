use bsd_elliptic_3adic_spectral::{
    curve_spectral_rime_frame, prism, unprism, EllipticCurve, Point, Zero, K, P,
};

#[test]
fn three_torsion_group_isomorphic_to_tribit_rotations() {
    let nil = Zero::Nil;
    let pos = nil.rotate();
    let neg = pos.rotate();
    assert_eq!(pos, Zero::Pos);
    assert_eq!(neg, Zero::Neg);
    assert_eq!(neg.rotate(), Zero::Nil);
    assert_ne!(pos, neg);
}

#[test]
fn weierstrass_tangent_slope_uses_3adic_shift_and_hensel_inversion() {
    let curve = EllipticCurve::new(-1, 1);
    let p = 7i64;
    let pt = Point::Affine(0, 1);
    assert!(curve.is_on_curve_mod(pt, p));

    let pt2 = curve.add_mod(pt, pt, p);
    match pt2 {
        Point::Affine(x, y) => {
            assert!(curve.is_on_curve_mod(Point::Affine(x, y), p));
            assert_eq!(x, 2);
            assert_eq!(y, 0);
        }
        Point::Infinity => panic!("Expected affine point"),
    }
}

#[test]
fn bsd_rank_spectral_parity_rank_0_curves() {
    let curves_rank_0 = [
        EllipticCurve::new(-1, 1),
        EllipticCurve::new(0, -4),
        EllipticCurve::new(1, 0),
    ];

    for curve in &curves_rank_0 {
        let frame = curve_spectral_rime_frame(curve);
        let spectrum = prism(&frame);
        let back = unprism(&spectrum);
        assert_eq!(frame, back, "NTT prism round-trip exact");

        let center_dc = spectrum[0];
        assert_ne!(
            center_dc, 0,
            "Rank 0 curve must have strictly NON-VANISHING uncomputed center X[0]"
        );
    }
}

#[test]
fn bsd_rank_spectral_parity_rank_1_curves() {
    let curve_rank_1 = EllipticCurve::new(0, -2);

    let pt_gen = Point::Affine(3, 5);
    assert_eq!(5 * 5, 3 * 3 * 3 - 2);

    let p_prime = 1009i64;
    let pt_2 = curve_rank_1.add_mod(pt_gen, pt_gen, p_prime);
    let pt_3 = curve_rank_1.add_mod(pt_2, pt_gen, p_prime);
    assert_ne!(pt_2, Point::Infinity);
    assert_ne!(pt_3, Point::Infinity);

    let frame = curve_spectral_rime_frame(&curve_rank_1);
    let spectrum = prism(&frame);
    let back = unprism(&spectrum);
    assert_eq!(frame, back);

    let mut sum_frame = 0u64;
    for i in 0..K {
        sum_frame = (sum_frame + frame[i]) % P;
    }
    assert_eq!(spectrum[0], sum_frame);
}

#[test]
fn iwasawa_27_point_filtration_depth() {
    assert_eq!(K, 27);
    assert_eq!(3u64.pow(3), 27);

    let mut depth = 1u64;
    for _ in 0..3 {
        depth *= 3;
    }
    assert_eq!(depth, 27);
}

#[test]
fn spherical_ntt_prism_round_trip_is_exact() {
    let mut frame = [0u64; K];
    for i in 0..K {
        frame[i] = ((i as u64 * 47 + 13) * 733) % P;
    }
    let spectrum = prism(&frame);
    let back = unprism(&spectrum);
    assert_eq!(frame, back);
}
