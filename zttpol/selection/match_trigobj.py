import law
from columnflow.selection import SelectionResult
from columnflow.util import maybe_import

from zttpol.util import filter_by_triggers, get_objs_p4, trigger_matching_extra, trigger_object_matching_deep

np = maybe_import("numpy")
ak = maybe_import("awkward")



def match_trigobjs_semilep(
        leps_pair: ak.Array,
        trigger_results: SelectionResult,
        **kwargs,
) -> tuple[ak.Array, ak.Array]:

    ch = kwargs.get("channel_1") # e,mu
    
    # extract the trigger names, types & others from trigger_results.x (aux)
    trigger_ids           = trigger_results.x.trigger_ids
    trigger_types         = trigger_results.x.trigger_types
    leg1_minpt            = trigger_results.x.leg1_minpt
    leg2_minpt            = trigger_results.x.leg2_minpt
    leg1_maxeta           = trigger_results.x.leg1_maxeta
    leg2_maxeta           = trigger_results.x.leg2_maxeta
    leg1_matched_trigobjs = trigger_results.x.leg1_matched_trigobjs
    leg2_matched_trigobjs = trigger_results.x.leg2_matched_trigobjs

    has_single_e_triggers = trigger_types == f"single_{ch}"
    has_cross_e_triggers  = trigger_types == f"cross_{ch}_tau"
    has_e_triggers = (has_single_e_triggers | has_cross_e_triggers)


    etau_pair  = filter_by_triggers(leps_pair, has_e_triggers)
    eles, taus = ak.unzip(etau_pair)

    
    # Event level masks
    # if events have etau
    has_etau_pairs = ak.fill_none(ak.num(eles, axis=1) > 0, False)

    mask_has_single_e_triggers_and_has_etau_pairs = has_single_e_triggers & has_etau_pairs
    mask_has_cross_e_triggers_and_has_etau_pairs  = has_cross_e_triggers & has_etau_pairs
    mask_has_e_triggers_and_has_etau_pairs = has_e_triggers & has_etau_pairs

    # filtering out the info based on the masks defined just above
    # for etau and mutau type of events, separate masks are created
    # for single and cross triggered events
    # for single e triggers
    single_etau_trigger_types          = trigger_types[mask_has_single_e_triggers_and_has_etau_pairs]
    single_etau_trigger_ids            = trigger_ids[mask_has_single_e_triggers_and_has_etau_pairs]
    single_etau_leg_1_minpt            = leg1_minpt[mask_has_single_e_triggers_and_has_etau_pairs]
    single_etau_leg_1_maxeta           = leg1_maxeta[mask_has_single_e_triggers_and_has_etau_pairs] 
    single_etau_leg_1_matched_trigobjs = leg1_matched_trigobjs[mask_has_single_e_triggers_and_has_etau_pairs]
    # for cross e triggers
    cross_etau_trigger_types          = trigger_types[mask_has_cross_e_triggers_and_has_etau_pairs]
    cross_etau_trigger_ids            = trigger_ids[mask_has_cross_e_triggers_and_has_etau_pairs]
    cross_etau_leg_1_minpt            = leg1_minpt[mask_has_cross_e_triggers_and_has_etau_pairs]
    cross_etau_leg_2_minpt            = leg2_minpt[mask_has_cross_e_triggers_and_has_etau_pairs] 
    cross_etau_leg_1_maxeta           = leg1_maxeta[mask_has_cross_e_triggers_and_has_etau_pairs]
    cross_etau_leg_2_maxeta           = leg2_maxeta[mask_has_cross_e_triggers_and_has_etau_pairs] 
    cross_etau_leg_1_matched_trigobjs = leg1_matched_trigobjs[mask_has_cross_e_triggers_and_has_etau_pairs]
    cross_etau_leg_2_matched_trigobjs = leg2_matched_trigobjs[mask_has_cross_e_triggers_and_has_etau_pairs]
    # concatenating single and cross info, so that it does not remain depend on the trigger hierarchy
    etau_trigger_types          = ak.concatenate([single_etau_trigger_types, cross_etau_trigger_types], axis=-1)
    etau_trigger_ids            = ak.concatenate([single_etau_trigger_ids, cross_etau_trigger_ids], axis=-1)


    # to convert the masks to event level
    # e.g. events with etau pair and pass electron triggers
    mask_has_single_e_triggers_and_has_e_evt_level = ak.any(mask_has_single_e_triggers_and_has_etau_pairs, axis=1)
    mask_has_cross_e_triggers_and_has_e_evt_level  = ak.any(mask_has_cross_e_triggers_and_has_etau_pairs, axis=1)
    mask_has_e_triggers_and_has_e_evt_level        = ak.any(mask_has_e_triggers_and_has_etau_pairs, axis=1)
    
    # dummy bool array
    trigobj_matched_mask_dummy = ak.from_regular((trigger_ids > 0)[:,:0][:,None])

    single_el_trigobj_matched_mask = trigger_object_matching_deep(eles,
                                                                  single_etau_leg_1_matched_trigobjs,
                                                                  single_etau_leg_1_minpt,
                                                                  single_etau_leg_1_maxeta,
                                                                  True)

    #print("ff")
    cross_el_trigobj_matched_mask_leg1 = trigger_object_matching_deep(eles,
                                                                      cross_etau_leg_1_matched_trigobjs,
                                                                      cross_etau_leg_1_minpt,
                                                                      cross_etau_leg_1_maxeta,
                                                                      True)
    #print("gg")
    cross_el_trigobj_matched_mask_leg2 = trigger_object_matching_deep(taus,
                                                                      cross_etau_leg_2_matched_trigobjs,
                                                                      cross_etau_leg_2_minpt,
                                                                      cross_etau_leg_2_maxeta,
                                                                      True)

    cross_el_trigobj_matched_mask = (cross_el_trigobj_matched_mask_leg1 & cross_el_trigobj_matched_mask_leg2)

    
    #from IPython import embed; embed()

    single_el_trigobj_matched_mask_evt_level = ak.fill_none(ak.firsts(ak.any(single_el_trigobj_matched_mask, axis=1), axis=1), False)
    cross_el_trigobj_matched_mask_evt_level = ak.fill_none(ak.firsts(ak.any(cross_el_trigobj_matched_mask, axis=1), axis=1), False) 
    

    match_single = mask_has_single_e_triggers_and_has_e_evt_level & single_el_trigobj_matched_mask_evt_level
    match_cross  = mask_has_cross_e_triggers_and_has_e_evt_level & cross_el_trigobj_matched_mask_evt_level

    
    el_trigobj_matched_mask = ak.where(match_single, single_el_trigobj_matched_mask, cross_el_trigobj_matched_mask)
    el_trigobj_matched_mask = ak.fill_none(ak.firsts(el_trigobj_matched_mask, axis=-1), False)

    new_eles = eles[el_trigobj_matched_mask]
    new_taus = taus[el_trigobj_matched_mask]

    #trigIds = ak.where(match_single, single_etau_trigger_ids, cross_etau_trigger_ids)
    #trigTypes = ak.where(match_single, single_etau_trigger_types, cross_etau_trigger_types)

    trigTypes = ak.concatenate([single_etau_trigger_types, cross_etau_trigger_types], axis=1)
    trigIds = ak.concatenate([single_etau_trigger_ids, cross_etau_trigger_ids], axis=1)
    ids = ak.values_astype(trigIds, 'int64')



    #from IPython import embed; embed()
    
    leps_pair = ak.zip([new_eles, new_taus])
    
    #from IPython import embed; embed()

    #ids_dummy = ak.from_regular((trigger_ids > 0)[:,:0])
    #ids = ak.where(mask_has_tau_triggers_and_has_tau_pairs_evt_level, ids, ids_dummy)

    return leps_pair, ids, trigTypes




def match_trigobjs_fullhad(
        leps_pair: ak.Array,
        trigger_results: SelectionResult,
        jets: ak.Array,
        **kwargs,
) -> tuple[ak.Array, ak.Array]:

    # extract the trigger names, types & others from trigger_results.x (aux)
    trigger_ids           = trigger_results.x.trigger_ids
    trigger_types         = trigger_results.x.trigger_types
    leg1_minpt            = trigger_results.x.leg1_minpt
    leg2_minpt            = trigger_results.x.leg2_minpt
    leg3_minpt            = trigger_results.x.leg3_minpt # NEW
    leg1_maxeta           = trigger_results.x.leg1_maxeta
    leg2_maxeta           = trigger_results.x.leg2_maxeta
    leg3_maxeta           = trigger_results.x.leg3_maxeta # NEW
    leg1_matched_trigobjs = trigger_results.x.leg1_matched_trigobjs
    leg2_matched_trigobjs = trigger_results.x.leg2_matched_trigobjs
    leg3_matched_trigobjs = trigger_results.x.leg3_matched_trigobjs # NEW

    has_tau_triggers     = ((trigger_types == "cross_tau_tau") | (trigger_types == "cross_tau_tau_jet"))

    old_leps_pair = leps_pair

    #from IPython import embed; embed()
    # ----- WARNING : TIME CONSUMING ------ #
    leps_pair = filter_by_triggers(leps_pair, has_tau_triggers) 

    taus1, taus2 = ak.unzip(leps_pair)

    #from IPython import embed; embed()

    #jets = filter_by_triggers(jets, trigger_types == "cross_tau_tau_jet") # NEW
    #jets = events.Jet[jet_indices]
    jets = filter_by_triggers(jets, has_tau_triggers) # NEW

    # IMPORTANT : APPLY JET PT > 60 GeV BEFORE MATCHING WITH JET LEG
    jets = jets[jets.pt > 60.0]

    
    # Event level masks
    # if events have tau
    has_tau_pairs = ak.fill_none(ak.num(taus1, axis=1) > 0, False)

    # events must be fired by tau triggers and there is ta inside
    mask_has_tau_triggers_and_has_tau_pairs = has_tau_triggers & has_tau_pairs

    tautau_trigger_ids            = trigger_ids[mask_has_tau_triggers_and_has_tau_pairs]    
    tautau_trigger_types          = trigger_types[mask_has_tau_triggers_and_has_tau_pairs]
    tautau_leg_1_minpt            = leg1_minpt[mask_has_tau_triggers_and_has_tau_pairs]
    tautau_leg_2_minpt            = leg2_minpt[mask_has_tau_triggers_and_has_tau_pairs]
    tautau_leg_3_minpt            = leg3_minpt[mask_has_tau_triggers_and_has_tau_pairs] # NEW
    tautau_leg_1_maxeta           = leg1_maxeta[mask_has_tau_triggers_and_has_tau_pairs]
    tautau_leg_2_maxeta           = leg2_maxeta[mask_has_tau_triggers_and_has_tau_pairs]
    tautau_leg_3_maxeta           = leg3_maxeta[mask_has_tau_triggers_and_has_tau_pairs] # NEW
    tautau_leg_1_matched_trigobjs = leg1_matched_trigobjs[mask_has_tau_triggers_and_has_tau_pairs]
    tautau_leg_2_matched_trigobjs = leg2_matched_trigobjs[mask_has_tau_triggers_and_has_tau_pairs]
    tautau_leg_3_matched_trigobjs = leg3_matched_trigobjs[mask_has_tau_triggers_and_has_tau_pairs] # NEW

    #from IPython import embed; embed()
    mask_has_tau_triggers_and_has_tau_pairs_evt_level = ak.fill_none(ak.any(mask_has_tau_triggers_and_has_tau_pairs, axis=1), False)
    jets_dummy = jets[:,:0]
    #jets = jets[mask_has_tau_triggers_and_has_tau_pairs]
    jets = ak.where(mask_has_tau_triggers_and_has_tau_pairs_evt_level, jets, jets_dummy)
    

    # tau1            : [ [    t11,        t12      ], [  t11  ] ]
    # tau2            : [ [    t21,        t22      ], [  t21  ] ]

    # trigger_ids     : [ [    1120,       1121     ], [  1120 ] ]
    # leg1/2_minpt    : [ [     40 ,        40      ], [   40  ] ]
    # leg1/2_trigobjs : [ [[o1,o2,o3], [o1,o2,o3,o4]], [[o1,o2]] ]

    # In [67]: taus1.pt[1044]
    # Out[67]: <Array [111, 111, 58.7] type='3 * float32'>

    # In [72]: ak.to_list(tautau_leg_1_matched_trigobjs.pt[1044])
    # Out[72]: 
    # [[110.203125, 53.0, 61.3828125, 210.625],
    #  [110.203125, 53.0, 61.3828125, 210.625, 35.0703125]]
    
    # In [70]: ak.to_list(ak.any(dr[1044] < 0.5, axis=-1))
    # Out[70]: [[True, True], [True, True], [True, True]]

    # same for tau1-legs2 and for tau2
    
    # final mask may look like this:
    # [[True, True], [True, True], [False, True]]

    # is tau1 matched to leg1?
    pass_tau1_leg1_triglevel = trigger_object_matching_deep(taus1,
                                                            tautau_leg_1_matched_trigobjs,
                                                            tautau_leg_1_minpt,
                                                            tautau_leg_1_maxeta,
                                                            True)
    # is tau1 matched to leg2?
    pass_tau1_leg2_triglevel = trigger_object_matching_deep(taus1,
                                                            tautau_leg_2_matched_trigobjs,
                                                            tautau_leg_2_minpt,
                                                            tautau_leg_2_maxeta,
                                                            True)

    # is tau2 matched to leg1?
    pass_tau2_leg1_triglevel = trigger_object_matching_deep(taus2,
                                                            tautau_leg_1_matched_trigobjs,
                                                            tautau_leg_1_minpt,
                                                            tautau_leg_1_maxeta,
                                                            True)
    # is tau2 matched to leg2?
    pass_tau2_leg2_triglevel = trigger_object_matching_deep(taus2,
                                                            tautau_leg_2_matched_trigobjs,
                                                            tautau_leg_2_minpt,
                                                            tautau_leg_2_maxeta,
                                                            True)

    # is any jet matched to leg3?
    pass_jets_leg3_triglevel, jets = trigger_object_matching_jets_deep(jets,
                                                                       tautau_leg_3_matched_trigobjs,
                                                                       tautau_leg_3_minpt,
                                                                       tautau_leg_3_maxeta,
                                                                       True)
    
    #from IPython import embed; embed()





    

    # tau1 to leg1 & tau2 to leg2 or, tau1 to leg2 & tau2 to leg1
    pass_taus_legs = (pass_tau1_leg1_triglevel & pass_tau2_leg2_triglevel) | (pass_tau1_leg2_triglevel & pass_tau2_leg1_triglevel)

    #pass_jet_leg_jet_level = ak.any(pass_jets_leg3_triglevel, axis=-1)
    #pass_jet_leg_jet_level_1 = pass_jet_leg_jet_level[:,:,0:1]
    #pass_jet_leg_jet_level_2 = pass_jet_leg_jet_level[:,:,1:2]
    #pass_jet_leg_any_1 = ak.any(pass_jet_leg_jet_level_1, axis=1)[:,None]
    #pass_jet_leg_any_2 = ak.any(pass_jet_leg_jet_level_2, axis=1)[:,None]
    #pass_jet_leg_tau_pair_level = ak.concatenate([pass_jet_leg_any_1, pass_jet_leg_any_2], axis=1)

    pass_jet_leg_jet_level = ak.any(pass_jets_leg3_triglevel, axis=-1)


    #dummy_true = ak.values_astype(ak.ones_like(tautau_trigger_ids), np.bool)
    #ak.where(tautau_trigger_ids == 15152, pass_jet_leg_trg_level, dummy_true)


    

    #pass_jet_leg_2 = ak.any(pass_jets_leg3_triglevel, axis=1)
    #temp = ak.values_astype(ak.ones_like(tautau_trigger_ids), np.bool)
    #dummy = temp[:,:0]
    #num_mask = ak.num(tautau_trigger_ids, axis=1) > 0
    #temp2 = ak.where(num_mask, temp[:,None], dummy)
    #pass_jet_leg_2 = ak.where(num_mask, pass_jet_leg[:,None], dummy)

    #from IPython import embed; embed()

    #pass_jet_leg2 = ak.where(num_mask, pass_jet_leg_2[:,None], dummy)
    #dummy_2 = tautau_trigger_ids[:,:0]
    #ids_2 = ak.where(num_mask, tautau_trigger_ids[:,None], dummy_2)
    
    
    #ditaujet_mask = (tautau_trigger_ids == 15152)
    #ak.where(ditaujet_mask, pass_jet_leg, temp2)
    
    
    #mask_has_tau_triggers_and_has_tau_pairs_evt_level = ak.fill_none(ak.any(mask_has_tau_triggers_and_has_tau_pairs, axis=1), False)
    trigobj_matched_mask_dummy = ak.from_regular((trigger_ids > 0)[:,:0][:,None])
    
    pass_taus_legs = ak.where(mask_has_tau_triggers_and_has_tau_pairs_evt_level, pass_taus_legs, trigobj_matched_mask_dummy)
    pass_taus_legs = ak.enforce_type(ak.values_astype(pass_taus_legs, "bool"), "var * var * bool") # 100000 * var * var * bool
    
    
    pass_taus = ak.fill_none(ak.any(pass_taus_legs, axis=-1), False)

    #tautau_trigger_ids_dummy = tautau_trigger_ids[:,:0]
    #tautau_trigger_ids = ak.where(mask_has_tau_triggers_and_has_tau_pairs_evt_level, tautau_trigger_ids, tautau_trigger_ids_dummy)
    
    
    tautau_trigger_ids_brdcst, _ = ak.broadcast_arrays(tautau_trigger_ids[:,None], pass_taus)
    #tautau_trigger_ids_brdcst = ak.where(mask_has_tau_triggers_and_has_tau_pairs_evt_level,
    #                                     tautau_trigger_ids_brdcst, 
    #                                     trigobj_matched_mask_dummy)

    tautau_trigger_ids_brdcst = tautau_trigger_ids_brdcst[pass_taus]
    # -------------> Order of triggers is important. Mind the hierarchy
    # this ids are the trig-obj match ids

    tautau_trigger_ids_brdcst_dummy = ak.from_regular(tautau_trigger_ids_brdcst[:,:0][:,None])
    #_tautau_trigger_ids_brdcst = ak.fill_none(ak.firsts(tautau_trigger_ids_brdcst, axis=1), 0)
    tautau_trigger_ids_brdcst = ak.where(ak.num(tautau_trigger_ids_brdcst) > 0, tautau_trigger_ids_brdcst, tautau_trigger_ids_brdcst_dummy)
    
    #ids = ak.fill_none(ak.firsts(tautau_trigger_ids_brdcst, axis=-1), -1)

    ids = ak.Array(ak.to_list(ak.firsts(tautau_trigger_ids_brdcst, axis=1))) # BAD Practice !!!
    
    #ids_dummy = ak.from_regular((trigger_ids > 0)[:,:0])
    #ids = ak.where(mask_has_tau_triggers_and_has_tau_pairs_evt_level, ids, ids_dummy)
    
    ids = ids[ak.fill_none(ak.firsts(pass_taus_legs, axis=1), False)]
    ids = ak.values_astype(ids, 'int64')

    new_taus1 = taus1[pass_taus]
    new_taus2 = taus2[pass_taus]
    #ids = ids[pass_taus]
    
    leps_pair = ak.zip([new_taus1, new_taus2])

    tautau_trigger_types_brdcst, _ = ak.broadcast_arrays(tautau_trigger_types[:,None], pass_taus)
    #tautau_trigger_ids_brdcst = ak.where(mask_has_tau_triggers_and_has_tau_pairs_evt_level,
    #                                     tautau_trigger_ids_brdcst, 
    #                                     trigobj_matched_mask_dummy)

    tautau_trigger_types_brdcst = tautau_trigger_types_brdcst[pass_taus]
    # -------------> Order of triggers is important. Mind the hierarchy
    # this ids are the trig-obj match ids


    #from IPython import embed; embed()

    tautau_trigger_types_brdcst_dummy = ak.from_regular(tautau_trigger_types_brdcst[:,:0][:,None])
    tautau_trigger_types_brdcst = ak.where(ak.num(tautau_trigger_types_brdcst) > 0, tautau_trigger_types_brdcst, tautau_trigger_types_brdcst_dummy)
    
    
    types = ak.Array(ak.to_list(ak.firsts(tautau_trigger_types_brdcst, axis=1))) # BAD Practice !!!
    types = types[ak.fill_none(ak.firsts(pass_taus_legs, axis=1), False)]
    
    #types = ak.fill_none(ak.firsts(tautau_trigger_types_brdcst, axis=-1), "")

    #from IPython import embed; embed()


    #pass_jet_leg_2 = pass_jet_leg[:,None]
    #dummy = pass_jet_leg_jet_level[:,:0]
    #pass_jet_leg_2 = ak.where(mask_has_tau_triggers_and_has_tau_pairs_evt_level, pass_jet_leg_2, dummy)

    #ids_mask_ditaujet = ak.sort(ids == 15152, ascending=False)
    #mask_ditaujet = ak.fill_none(ak.firsts(ids_mask_ditaujet, axis=1), False)
    #ids_mask_ditau = ids == 15151
    #mask_ditau = ak.fill_none(ak.firsts(ids_mask_ditaujet, axis=1), False)
    #mask_ditaujet_with_jetmatch = mask_ditaujet & pass_jet_leg

    #mask_all_matched = mask_ditau | mask_ditaujet_with_jetmatch

    #ids_ditau = 

    #from IPython import embed; embed()


    
    #-----ids = ak.enforce_type(ids, "var * int64")
    #-----types = ak.enforce_type(types, "var * string")
    """
    jets_idx = jets[pass_jet_leg_jet_level].rawIdx[:,:1]
    """
    jets_idx = jets.rawIdx[:,:1]
    #jets_idx = ak.local_index(jets.pt)[pass_jet_leg_jet_level][:,:1]
    
    return leps_pair, ids, types, jets_idx
