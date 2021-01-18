
#' Returns a list of 39 tibbles, with the name and data returned from each respective reference function.
#' @return Returns a list of 39 tibbles, with the name and data returned from each respective reference function.
#' @param tidy Fix the variable names in each tibble to remove special characters and superfluous text, and convert all the variable names to snake_case. Defaults to TRUE.
#' @param tidy_style The style to convert variable names to, if tidy=TRUE. Accepts one of "snake_case", "camelCase" and "period.case". Defaults to "snake_case".
#' @seealso \code{\link{mnis_reference}}
#' @rdname mnis_all_reference
#' @export
#' @examples \dontrun{
#'
#' z <- mnis_all_reference()
#'
#'}

mnis_all_reference <- function(tidy = TRUE, tidy_style="snake_case") {

    address_types <- ref_address_types(tidy = tidy, tidy_style="snake_case")

    answering_bodies <- ref_answering_bodies(tidy = tidy, tidy_style="snake_case")

    area_types <- ref_area_types(tidy = tidy, tidy_style="snake_case")

    areas <- ref_areas(tidy = tidy, tidy_style="snake_case")

    biography_categories <- ref_biography_categories(tidy = tidy, tidy_style="snake_case")

    cabinets <- ref_cabinets(tidy = tidy, tidy_style="snake_case")

    committee_types <- ref_committee_types(tidy = tidy, tidy_style="snake_case")

    committees <- ref_committees(tidy = tidy, tidy_style="snake_case")

    constituencies <- ref_constituencies(tidy = tidy, tidy_style="snake_case")

    constituency_areas <- ref_constituency_areas(tidy = tidy, tidy_style="snake_case")

    constituency_types <- ref_constituency_types(tidy = tidy, tidy_style="snake_case")

    countries <- ref_countries(tidy = tidy, tidy_style="snake_case")

    departments <- ref_departments(tidy = tidy, tidy_style="snake_case")

    disqualification_types <- ref_disqualification_types(tidy = tidy, tidy_style="snake_case")

    election_types <- ref_election_types(tidy = tidy, tidy_style="snake_case")

    elections <- ref_elections(tidy = tidy, tidy_style="snake_case")

    end_reasons <- ref_end_reasons(tidy = tidy, tidy_style="snake_case")

    experience_types <- ref_experience_types(tidy = tidy, tidy_style="snake_case")

    government_post_departments <- ref_government_post_departments(tidy = tidy, tidy_style="snake_case")

    government_posts <- ref_government_posts(tidy = tidy, tidy_style="snake_case")

    government_ranks <- ref_government_ranks(tidy = tidy, tidy_style="snake_case")

    honour_lists <- ref_honour_lists(tidy = tidy, tidy_style="snake_case")

    honourary_prefixes <- ref_honourary_prefixes(tidy = tidy, tidy_style="snake_case")

    honours <- ref_honours(tidy = tidy, tidy_style="snake_case")

    interest_categories <- ref_interest_categories(tidy = tidy, tidy_style="snake_case")

    lords_membership_types <- ref_lords_membership_types(tidy = tidy, tidy_style="snake_case")

    lords_ranks <- ref_lords_ranks(tidy = tidy, tidy_style="snake_case")

    opposition_post_departments <- ref_opposition_post_departments(tidy = tidy, tidy_style="snake_case")

    opposition_posts <- ref_opposition_posts(tidy = tidy, tidy_style="snake_case")

    opposition_ranks <- ref_opposition_ranks(tidy = tidy, tidy_style="snake_case")

    other_parliaments <- ref_other_parliaments(tidy = tidy, tidy_style="snake_case")

    parliament_types <- ref_parliament_types(tidy = tidy, tidy_style="snake_case")

    parliamentary_posts <- ref_parliamentary_posts(tidy = tidy, tidy_style="snake_case")

    parliamentary_ranks <- ref_parliamentary_ranks(tidy = tidy, tidy_style="snake_case")

    parties <- ref_parties(tidy = tidy, tidy_style="snake_case")

    party_sub_types <- ref_party_sub_types(tidy = tidy, tidy_style="snake_case")

    photo_outputs <- ref_photo_outputs(tidy = tidy, tidy_style="snake_case")

    statuses <- ref_statuses(tidy = tidy, tidy_style="snake_case")

    titles <- ref_titles(tidy = tidy, tidy_style="snake_case")

    ## make the reference

    ref_list <- list()

    ref_list <- list(address_types = address_types,
                     answering_bodies = answering_bodies,
                     area_types = area_types,
                     areas = areas,
                     biography_categories = biography_categories,
                     cabinets = cabinets,
                     committee_types = committee_types,
                     committees = committees,
                     constituencies = constituencies,
                     constituency_areas = constituency_areas,
                     constituency_types = constituency_types,
                     countries = countries,
                     departments = departments,
                     disqualification_types = disqualification_types,
                     election_types = election_types,
                     elections = elections,
                     end_reasons = end_reasons,
                     experience_types = experience_types,
                     government_post_departments = government_post_departments,
                     government_posts = government_posts,
                     government_ranks = government_ranks,
                     honour_lists = honour_lists, honourary_prefixes = honourary_prefixes,
                     honours = honours,
                     interest_categories = interest_categories,
                     lords_membership_types = lords_membership_types,
                     lords_ranks = lords_ranks,
                     opposition_post_departments = opposition_post_departments,
                     opposition_posts = opposition_posts,
                     opposition_ranks = opposition_ranks,
                     other_parliaments = other_parliaments,
                     parliament_types = parliament_types,
                     parliamentary_posts = parliamentary_posts,
                     parliamentary_ranks = parliamentary_ranks,
                     parties = parties,
                     party_sub_types = party_sub_types,
                     photo_outputs = photo_outputs,
                     statuses = statuses,
                     titles = titles)

    ref_list

}
