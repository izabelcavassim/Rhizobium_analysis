import itertools
from copy import deepcopy, copy
from operator import itemgetter

import numpy as np

__author__ = 'Asger Bachmann (agb@birc.au.dk)'


class Gene:
    def __init__(self, name, seq_name, start, end, reverse):
        self.name = name
        self.seq_name = seq_name
        self.start = start
        self.end = end
        self.reverse = reverse
        self.neighbor_l = None
        self.neighbor_r = None


def get_gene_distance(g1, g2):
    if g1.seq_name != g2.seq_name:
        return None
    elif g1.start < g2.start:
        return g2.start - g1.end
    else:
        return g1.start - g2.end


def get_min_dist_to_neighbors(n_map, unambiguous_strain, unambiguous_gene, genes):
    n_map = n_map[unambiguous_strain] if unambiguous_strain in n_map else []
    min_dist = None
    for n_map_g in n_map:
        neighbour_gene = genes[unambiguous_strain][n_map_g]
        dist = get_gene_distance(unambiguous_gene, neighbour_gene)
        if dist is not None:
            min_dist = min(min_dist, dist) if min_dist is not None else dist
    return min_dist


def all_unique(x):
    seen = set()
    return not any(i in seen or seen.add(i) for i in x)


def get_groups(base_groups):
    groups = []
    unambiguous_groups = []
    ambiguous_groups = []
    member_mapping = {}
    for group_id, group in base_groups:
        unambiguous = True
        group_members = {}
        for member_strain, member_gene in group:
            if member_strain in group_members:
                unambiguous = False
                group_members[member_strain].append(member_gene)
            else:
                group_members[member_strain] = [member_gene]

        for s, gs in group_members.items():
            if s not in member_mapping:
                member_mapping[s] = {}
            member_mapping[s].update({g: group_id for g in gs})  # group_members for g in gs})

        groups.append((group_id, group_members))
        if unambiguous:
            unambiguous_groups.append((group_id, group_members))
        else:
            ambiguous_groups.append((group_id, group_members))

    return unambiguous_groups, ambiguous_groups, groups, member_mapping


def get_disambiguated_groups2(base_groups, genes, neighborhood_size=2, top_matchers_threshold=0.0,
                              gene_inclusion_threshold=0.0, verbose=False):
    unambiguous_groups, ambiguous_groups, groups, member_mapping = get_groups(base_groups)

    if verbose:
        print("Number of groups prior to process: {}".format(len(groups)))

    resulting_groups = []

    for group_name, members in groups:
        print("Processing {}...".format(group_name))
        # Ensure reproducibility by sorting the members to keep a stable iteration order.
        flattened_members = sorted([(s, g) for s, gs in members.items() for g in gs])
        score_matrix = np.zeros((len(flattened_members), len(flattened_members)))
        individual_scores = np.zeros(len(flattened_members))
        genes_to_skip = set()

        for i, (s1, g1) in enumerate(flattened_members):
            g1 = genes[s1][g1]
            ns1 = get_neighbors(g1, s1, member_mapping, neighborhood_size)

            for j, (s2, g2) in enumerate(flattened_members[(i + 1):]):
                g2 = genes[s2][g2]
                ns2 = get_neighbors(g2, s2, member_mapping, neighborhood_size)

                overlap = ns1.intersection(ns2)
                score = len(overlap) / (neighborhood_size * 2) if len(overlap) > 0 else 0
                score_matrix[i][i + 1 + j] = score
                score_matrix[i + 1 + j][i] = score
                individual_scores[i] += score
                individual_scores[i + 1 + j] += score

        if max(individual_scores) != 0:
            individual_scores = individual_scores / max(individual_scores)
            genes_to_skip.update(np.where(individual_scores < gene_inclusion_threshold)[0])
        else:
            genes_to_skip.update(range(len(individual_scores)))

        # For each individual i, select the other individuals j with which i matches best.
        top_matchers = [set(j for j in range(len(score_matrix))
                            if j not in genes_to_skip and (score_matrix[i][j] > top_matchers_threshold or i == j))
                        if i not in genes_to_skip else set()
                        for i in range(len(score_matrix))]

        for i in range(len(score_matrix)):
            if len(top_matchers[i]) == 0:
                continue
            for j in range(i + 1, len(score_matrix)):
                if len(top_matchers[j]) == 0:
                    continue

                u = top_matchers[i].union(top_matchers[j])
                if len(top_matchers[i]) + len(top_matchers[j]) != len(u):
                    top_matchers[i] = u
                    top_matchers[j] = u

        separate_groups = set(frozenset(s) for s in top_matchers if len(s) > 1)
        for i, separate_group in enumerate(separate_groups):
            new_group_members = {}
            for m in separate_group:
                s, g = flattened_members[m]
                if s not in new_group_members:
                    new_group_members[s] = []
                new_group_members[s].append(g)

            resulting_groups.append(("{}_{}".format(group_name, i), new_group_members))

    if verbose:
        print("Number of groups after process: {}".format(len(resulting_groups)))

    return resulting_groups


def get_neighbors(g1, s1, member_mapping, neighborhood_size):
    # Follow the neighbor chain for up to neighborhood_size steps.
    result = set()
    n = g1.neighbor_l
    for i in range(neighborhood_size):
        if n is None:
            break

        if n.name in member_mapping[s1]:
            n_map = member_mapping[s1][n.name]
            if n_map is not None:
                result.add(n_map)

        n = n.neighbor_l

    n = g1.neighbor_r
    for i in range(neighborhood_size):
        if n is None:
            break

        if n.name in member_mapping[s1]:
            n_map = member_mapping[s1][n.name]
            if n_map is not None:
                result.add(n_map)

        n = n.neighbor_r
    return result


def get_disambiguated_groups(base_groups, genes, neighborhood_size=2, verbose=False):
    # The ratio of two most chosen path's counts.
    param_unique_path_count_ratio = 1.5
    # The ratio of two most chosen path's size.
    param_unique_path_size_ratio = 0.5
    # The ratio of two most chosen gene's weights.
    param_gene_weight_ratio = 1.5
    # The size of an ambiguous strain that we allow to split on.
    param_split_size = 2

    unambiguous_groups, ambiguous_groups, groups, member_mapping = get_groups(base_groups)

    while True:
        members_changed = False

        if verbose:
            print("Num groups: {}".format(len(groups)))
            print("\tNum unambiguous groups: {}".format(len(unambiguous_groups)))
            print("\tNum ambiguous groups: {}".format(len(ambiguous_groups)))

        groups_to_change = []
        groups_to_remove = []
        for g_idx, (group_name, members) in zip(range(len(ambiguous_groups)), ambiguous_groups):
            print("Processing group {}...".format(group_name))
            unambiguous_members = [(s, genes[s][gs[0]]) for s, gs in members.items() if len(gs) == 1]
            ambiguous_members = {s: gs for s, gs in members.items() if len(gs) > 1}

            num_unambiguous_members = len(unambiguous_members)
            # If a strain X has param_split_size members, split the group into param_split_size groups, each with 1
            # member from X. Choose the best one according to some measure of group member stability. Do this for all
            # strains with param_split_size members.
            split_candidates = [(s, gs) for s, gs in ambiguous_members.items()
                                if 1 < len(gs) < (param_split_size + 1)]

            if 0 < len(split_candidates):
                group_candidates = []
                for split_group_i, (split_candidate_strain, split_candidate_genes) in enumerate(split_candidates):
                    # print("Splitting {} in {} on strain {}...".format(group_name, len(split_candidate_genes),
                    #                                                   split_candidate_strain)
                    for split_candidate_gene in split_candidate_genes:
                        group_ambiguous_members = deepcopy(ambiguous_members)
                        del group_ambiguous_members[split_candidate_strain]
                        group_unambiguous_members = copy(unambiguous_members)
                        group_unambiguous_members.append(
                                (split_candidate_strain, genes[split_candidate_strain][split_candidate_gene]))
                        group_candidates.append((split_group_i, group_unambiguous_members, group_ambiguous_members))
            elif num_unambiguous_members == 0:
                # print("Ignoring group {} as it has no unambiguous strains and no split candidates.".format(group_name)
                continue
            else:
                # There is no need to split this group.
                group_candidates = [(None, unambiguous_members, ambiguous_members)]

            # Monitor the 'paths' chosen down through the ambiguous strains. If two or more separate paths exists with
            # all ambiguous strains in them, we can safely split the group into more groups per those paths.
            # For example, if we have a group, {A1,A2,B1,B2,C1,C2} with two paths, A1->B1->C1 and A2->B2->C2, we can
            # split the group into the two groups {A1,B1,C1} and {A2,B2,C2}.
            all_paths_chosen = {}
            # Similarly, count the number of times a gene is chosen for an ambiguous strain. The genes chosen the most
            # can be assumed to be the 'most correct' gene among the genes on the ambiguous strain.
            # For example, if 4 groups has A1 from strain S and only 1 group has A2 from strain S, we use A1.
            all_groups_chosen_gene_count = {}
            all_groups_chosen_gene_weight = {}
            for group_suffix, unambiguous_members, ambiguous_members in group_candidates:

                # Using each unambiguous strain as a reference strain, find the most stable group, i.e. the group
                # members that most reference strains agrees on.
                for reference_unambiguous_strain, reference_unambiguous_gene in unambiguous_members:
                    # Start the path with the unambiguous members.
                    path_chosen = [(unambiguous_strain, unambiguous_gene.name) for unambiguous_strain, unambiguous_gene
                                   in unambiguous_members]
                    total_min_dist = 0

                    for ambiguous_strain, ambiguous_genes in ambiguous_members.items():
                        min_dist = None
                        min_dist_gene = None
                        min_dist_gene_weight = 0

                        for g in ambiguous_genes:
                            ambiguous_gene = genes[ambiguous_strain][g]
                            l_dists = []
                            r_dists = []

                            # Follow the neighbor chain for up to neighborhood_size steps.
                            n_l = ambiguous_gene.neighbor_l
                            for n_l_idx in range(neighborhood_size):
                                if n_l is None:
                                    break

                                if n_l.name in member_mapping[ambiguous_strain]:
                                    n_l_map = member_mapping[ambiguous_strain][n_l.name]
                                    if n_l_map is not None:
                                        l_dist_part = get_min_dist_to_neighbors(n_l_map, reference_unambiguous_strain,
                                                                                reference_unambiguous_gene, genes)
                                        if l_dist_part is not None:
                                            l_dists.append(l_dist_part)

                                n_l = n_l.neighbor_l

                            n_r = ambiguous_gene.neighbor_r
                            for n_r_idx in range(neighborhood_size):
                                if n_r is None:
                                    break

                                if n_r.name in member_mapping[ambiguous_strain]:
                                    n_r_map = member_mapping[ambiguous_strain][n_r.name]
                                    if n_r_map is not None:
                                        r_dist_part = get_min_dist_to_neighbors(n_r_map, reference_unambiguous_strain,
                                                                                reference_unambiguous_gene, genes)
                                        if r_dist_part is not None:
                                            r_dists.append(r_dist_part)

                                n_r = n_r.neighbor_r

                            l_dist = None if len(l_dists) == 0 else min(l_dists)
                            r_dist = None if len(r_dists) == 0 else min(r_dists)

                            if l_dist is not None and r_dist is not None:
                                total_dist = (l_dist + r_dist) / 2
                            elif l_dist is not None:
                                total_dist = l_dist
                            elif r_dist is not None:
                                total_dist = r_dist
                            else:
                                total_dist = None

                            if total_dist is not None and (min_dist is None or total_dist < min_dist):
                                min_dist = total_dist
                                min_dist_gene = g
                                # We weight each gene by the number of neighbors that are actually correctly mapped. If
                                # two genes are chosen an equal number of times, we use this weight to choose between
                                # them.
                                min_dist_gene_weight = len(l_dists) + len(r_dists)

                        if min_dist_gene is not None:
                            total_min_dist += min_dist

                            if ambiguous_strain not in all_groups_chosen_gene_count:
                                all_groups_chosen_gene_count[ambiguous_strain] = {}
                                all_groups_chosen_gene_weight[ambiguous_strain] = {}
                            if min_dist_gene not in all_groups_chosen_gene_count[ambiguous_strain]:
                                all_groups_chosen_gene_count[ambiguous_strain][min_dist_gene] = 1
                                all_groups_chosen_gene_weight[ambiguous_strain][min_dist_gene] = min_dist_gene_weight
                            else:
                                all_groups_chosen_gene_count[ambiguous_strain][min_dist_gene] += 1
                                all_groups_chosen_gene_weight[ambiguous_strain][min_dist_gene] += min_dist_gene_weight

                            path_chosen.append((ambiguous_strain, min_dist_gene))

                    path_chosen = frozenset(path_chosen)
                    if path_chosen in all_paths_chosen:
                        all_paths_chosen[path_chosen] += 1
                    else:
                        all_paths_chosen[path_chosen] = 1

            # If the most chosen paths are non-overlapping, create new groups from these paths.
            if len(all_paths_chosen) > 1:
                most_chosen_path_count = max(all_paths_chosen.items(), key=itemgetter(1))[1]
                paths_with_matching_counts = [p for p, c in all_paths_chosen.items() if c == most_chosen_path_count]

                if len(paths_with_matching_counts) > 1:
                    if not any(len(p) < len(members) for p in paths_with_matching_counts):
                        flat_paths_with_matching_counts = itertools.chain.from_iterable(paths_with_matching_counts)
                        if all_unique(flat_paths_with_matching_counts):
                            if verbose:
                                print("Splitting group {} in {} non-overlapping paths...".format(
                                    group_name, len(paths_with_matching_counts)))

                            # Create new groups based on paths_with_matching_counts.
                            for p_id, p in zip(range(len(paths_with_matching_counts)), paths_with_matching_counts):
                                if len(p) == 1:
                                    if verbose:
                                        print("\tSplit of group resulted in a group of size 1. Removing the new group.")

                                    singleton_s, singleton_g = list(p)[0]
                                    if singleton_s in member_mapping and singleton_g in member_mapping[singleton_s]:
                                        del member_mapping[singleton_s][singleton_g]
                                    continue

                                new_group_members = {p_s: [p_g] for p_s, p_g in p}
                                new_group = ("{}_{}".format(group_name, p_id), new_group_members)
                                groups.append(new_group)
                                unambiguous_groups.append(new_group)

                                for new_s, new_gs in new_group_members.items():
                                    member_mapping[new_s].update({new_g: new_group_members for new_g in new_gs})
                            groups_to_remove.append(g_idx)
                            continue

            # If the most chosen path goes through all strains and is chosen more than the second-most chosen path
            # (according to a ratio), use the genes from that path.
            unique_path_found = False
            if len(all_paths_chosen) > 0:
                sorted_paths = sorted(all_paths_chosen.items(), key=itemgetter(1), reverse=True)
                most_chosen_path, most_chosen_path_count = sorted_paths[0]
                second_most_chosen_path, second_most_chosen_path_count = sorted_paths[1] if len(sorted_paths) > 1 \
                    else ([], 0)

                if len(most_chosen_path) == len(members):
                    if len(second_most_chosen_path) / len(most_chosen_path) < param_unique_path_size_ratio or \
                                    most_chosen_path_count > second_most_chosen_path_count * param_unique_path_count_ratio:
                        unique_path_found = True
                        for path_strain, path_gene in most_chosen_path:
                            members[path_strain] = [path_gene]
                            members_changed = True

            # Choose the most chosen genes as members of the groups.
            if not unique_path_found:
                discarded_genes_group_members = {}
                discarded_genes_group_unambiguous = True

                for gene_count_strain, gene_counts in all_groups_chosen_gene_count.items():
                    if len(gene_counts) == 0:
                        continue

                    sorted_gene_counts = sorted(gene_counts.items(), key=itemgetter(1), reverse=True)
                    best_gene_count = sorted_gene_counts[0][1]

                    if (num_unambiguous_members > 1 or len(group_candidates) > 1) and best_gene_count == 1:
                        continue

                    # Find all genes with the maximum count. If there are more than one, compare the genes by their
                    # weight.
                    best_gene_count_weights = [(n, all_groups_chosen_gene_weight[gene_count_strain][n]) for n, c in
                                               sorted_gene_counts if c == best_gene_count]
                    if len(best_gene_count_weights) > 1:
                        best_gene_count_weights = sorted(best_gene_count_weights, key=itemgetter(1), reverse=True)
                        # The most weighted gene must be weighted high enough compared to second most weighted gene.
                        if not best_gene_count_weights[0][1] >= param_gene_weight_ratio * best_gene_count_weights[1][1]:
                            continue

                    best_gene = best_gene_count_weights[0][0]
                    if len(gene_counts) > 1:
                        discarded_genes = [g_c[1] for g_c in sorted_gene_counts[1:]]
                        if len(discarded_genes) > 1:
                            discarded_genes_group_unambiguous = False
                        discarded_genes_group_members[gene_count_strain] = [discarded_genes]

                    # Update the group with the newly-found most plausible genes.
                    members[gene_count_strain] = [best_gene]
                    members_changed = True

                # Keep the discarded genes in a separate group. Their sequences are after all very similar.
                discarded_genes_group = ("{}_d".format(group_name), discarded_genes_group_members)
                groups.append(discarded_genes_group)
                if discarded_genes_group_unambiguous:
                    unambiguous_groups.append(discarded_genes_group)
                else:
                    ambiguous_groups.append(discarded_genes_group)

            ambiguous_members = [(s, gs) for s, gs in members.items() if len(gs) > 1]
            num_ambiguous_members = len(ambiguous_members)
            if num_ambiguous_members == 0:
                groups_to_change.append(g_idx)

        groups_changed = 0
        for g_idx in groups_to_change:
            g = ambiguous_groups.pop(g_idx - groups_changed)
            unambiguous_groups.append(g)
            groups_changed += 1

        # The ith entry of smaller_idx_groups_changed is the number of values in groups_to_change that are smaller than
        # the ith entry of groups_to_remove. We use these values to offset the index used to remove groups.
        smaller_idx_groups_changed = [0]
        current_remove_idx = 0
        for idx, g_idx in enumerate(groups_to_change):
            if current_remove_idx < len(groups_to_remove) and g_idx < groups_to_remove[current_remove_idx]:
                smaller_idx_groups_changed[current_remove_idx] += 1
            else:
                if current_remove_idx + 1 >= len(groups_to_remove):
                    break
                smaller_idx_groups_changed.append(smaller_idx_groups_changed[current_remove_idx] + 1)
                current_remove_idx += 1
        # Because groups_to_remove may be larger than groups_to_change, fill out the missing values.
        for idx in range(current_remove_idx, len(groups_to_remove)):
            smaller_idx_groups_changed.append(smaller_idx_groups_changed[current_remove_idx])

        groups_removed = 0
        for idx, g_idx in enumerate(groups_to_remove):
            g = ambiguous_groups.pop(g_idx - groups_removed - smaller_idx_groups_changed[idx])
            groups.remove(g)
            groups_removed += 1

        if not members_changed and groups_changed == 0 and groups_removed == 0:
            break

    return groups, ambiguous_groups, unambiguous_groups
