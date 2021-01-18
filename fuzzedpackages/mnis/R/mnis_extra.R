

#' A wrapper for \code{\link{mnis_additional}} functions.
#'
#' It combines the various options of mnis_additional into one dataframe, and the default is similar to \code{\link{mnis_full_biog}}. Variable descriptions are taken from the mnis website: <http://data.parliament.uk/membersdataplatform/memberquery.aspx>.
#'
#' @param ID The ID number of the member. Defaults to NULL. If NULL, returns
#' @param ref_dods Request based on the DODS membership ID scheme. Defaults to FALSE. If FALSE, requests data based on the default membership ID scheme.
#' @param addresses Member address information (e.g. website, twitter, consituency address etc...). Defaults to TRUE. If TRUE, address details are included in the tibble.
#' @param biography_entries Member biographical information (e.g. countries of interest, policy expertise etc...) Defaults to TRUE. If TRUE, biographical details are included in the tibble.
#' @param committees Committees a Member sits or has sat on as well details on committee chairing. Defaults to TRUE. If TRUE, committee details are included in the tibble.
#' @param constituencies constituencies a Member has represented. Defaults to TRUE. If TRUE, constituency details are included in the tibble.
#' @param elections_contested Elections a Member has contested but not won. Defaults to TRUE. If TRUE, details of unsuccessful election contests are included in the tibble.
#' @param experiences Non-parliamentary experience of a Member. Defaults to TRUE. If TRUE, extra-parliamentary experience details are included in the tibble.
#' @param government_posts Government posts a Member currently holds. Defaults to TRUE. If TRUE, government posts details are included in the tibble.
#' @param honours Honours (e.g. MBE, OBE etc...) held by a Member. Defaults to TRUE. If TRUE, honours details are included in the tibble.
#' @param house_memberships House membership list of a Member. Defaults to TRUE. If TRUE, house membership details are included in the tibble.
#' @param interests Registered interests (financial) of a Member. Defaults to TRUE. If TRUE, interest details are included in the tibble.
#' @param known_as Details of names a Member has chosen to be known as instead of their full title (House of Lords members only). Defaults to TRUE. If TRUE, known as details are included in the tibble.
#' @param maiden_speeches Maiden speech dates for a Member. Defaults to TRUE. If TRUE, maiden speech details are included in the tibble.
#' @param opposition_posts Opposition posts a Member has held. Defaults to TRUE. If TRUE, opposition post details are included in the tibble.
#' @param other_parliaments Other Parliaments that a Member has held a membership of. Defaults to TRUE. If TRUE, details of other parliaments are included in the tibble.
#' @param parliamentary_posts Parliamentary posts a Member has held. Defaults to TRUE. If TRUE, parliamentary posts details are included in the tibble.
#' @param parties Party affiliations of a Member. Defaults to TRUE. If TRUE, address details are included in the tibble.
#' @param preferred_names Full set of data about a Members' name (e.g. surname, forename, Honorary prefixes, full details of HoL title and rank etc...). Defaults to TRUE. If TRUE, preferred names details are included in the tibble.
#' @param staff The staff employed by a Member. Defaults to TRUE. If TRUE, staff details are included in the tibble.
#' @param statuses Status history (e.g. suspensions and disqualifications) for a Member. Defaults to TRUE. If TRUE, status details are included in the tibble.
#' @param tidy Fix the variable names in the tibble to remove special characters and superfluous text, and converts the variable names to a consistent style. Defaults to TRUE.
#' @param tidy_style The style to convert variable names to, if tidy=TRUE. Accepts one of "snake_case", "camelCase" and "period.case". Defaults to "snake_case"
#' @keywords mnis
#' @return A tibble with the requested data on a given MP.
#' @examples \dontrun{
#'
#' x <- mnis_extra(172)
#'
#' }
#' @export
#' @rdname mnis_extra
#' @seealso \code{\link{mnis_full_biog}} \code{\link{mnis_basic_details}} \code{\link{mnis_additional}}

mnis_extra <- function(ID, ref_dods = FALSE, addresses = TRUE, biography_entries = TRUE, committees = TRUE, constituencies = TRUE, elections_contested = TRUE, experiences = TRUE, government_posts = TRUE, honours = TRUE, house_memberships = TRUE, interests = TRUE, known_as = TRUE, maiden_speeches = TRUE, opposition_posts = TRUE, other_parliaments = TRUE, parliamentary_posts = TRUE, parties = TRUE, preferred_names = TRUE, staff = TRUE, statuses = TRUE, tidy = TRUE, tidy_style="snake_case") {

    ID <- as.character(ID)

    if(is.null(ID)==TRUE){
      stop("ID cannot be null", call. = FALSE)
    }

    mnis_df <- tibble::tibble(member_id=ID)

    if (addresses == TRUE) {
        addresses_df <- mnis_addresses(ID=ID, ref_dods = ref_dods, tidy=TRUE, tidy_style = "snake_case")
    suppressMessages(suppressMessages(suppressWarnings(mnis_df <- dplyr::inner_join(mnis_df, addresses_df))))}

    if (biography_entries == TRUE) {
        biography_entries_df <- mnis_biography_entries(ID=ID, ref_dods = ref_dods, tidy=TRUE, tidy_style = "snake_case")
    suppressMessages(suppressWarnings(mnis_df <- dplyr::inner_join(mnis_df, biography_entries_df)))}

    if (committees == TRUE) {
        committees_df <- mnis_committees(ID=ID, ref_dods = ref_dods, tidy=TRUE, tidy_style = "snake_case")
    suppressMessages(suppressWarnings(mnis_df <- dplyr::inner_join(mnis_df, committees_df)))}

    if (constituencies == TRUE) {
        constituencies_df <- mnis_constituencies(ID=ID, ref_dods = ref_dods, tidy=TRUE, tidy_style = "snake_case")
    suppressMessages(suppressWarnings(mnis_df <- dplyr::inner_join(mnis_df, constituencies_df)))}

    if (elections_contested == TRUE) {
        elections_contested_df <- mnis_elections_contested(ID=ID, ref_dods = ref_dods, tidy=TRUE, tidy_style = "snake_case")
    suppressMessages(suppressWarnings(mnis_df <- dplyr::inner_join(mnis_df, elections_contested_df)))}

    if (experiences == TRUE) {
        experiences_df <- mnis_experiences(ID=ID, ref_dods = ref_dods, tidy=TRUE, tidy_style = "snake_case")
    suppressMessages(suppressWarnings(mnis_df <- dplyr::inner_join(mnis_df, experiences_df)))}

    if (government_posts == TRUE) {
        government_posts_df <- mnis_government_posts(ID=ID, ref_dods = ref_dods, tidy=TRUE, tidy_style = "snake_case")
    suppressMessages(suppressWarnings(mnis_df <- dplyr::inner_join(mnis_df, government_posts_df)))}

    if (honours == TRUE) {
        honours_df <- mnis_honours(ID=ID, ref_dods = ref_dods, tidy=TRUE, tidy_style = "snake_case")
    suppressMessages(suppressWarnings(mnis_df <- dplyr::inner_join(mnis_df, honours_df)))}

    if (house_memberships == TRUE) {
        house_memberships_df <- mnis_house_memberships(ID=ID, ref_dods = ref_dods, tidy=TRUE, tidy_style = "snake_case")
    suppressMessages(suppressWarnings(mnis_df <- dplyr::inner_join(mnis_df, house_memberships_df)))}

    if (interests == TRUE) {
        interests_df <- mnis_interests(ID=ID, ref_dods = ref_dods, tidy=TRUE, tidy_style = "snake_case")
    suppressMessages(suppressWarnings(mnis_df <- dplyr::inner_join(mnis_df, interests_df)))}

    if (known_as == TRUE) {
        known_as_df <- mnis_known_as(ID=ID, ref_dods = ref_dods, tidy=TRUE, tidy_style = "snake_case")
    suppressMessages(suppressWarnings(mnis_df <- dplyr::inner_join(mnis_df, known_as_df)))}

    if (maiden_speeches == TRUE) {
        maiden_speeches_df <- mnis_maiden_speeches(ID=ID, ref_dods = ref_dods, tidy=TRUE, tidy_style = "snake_case")
    suppressMessages(suppressWarnings(mnis_df <- dplyr::inner_join(mnis_df, maiden_speeches_df)))}

    if (opposition_posts == TRUE) {
        opposition_posts_df <- mnis_opposition_posts(ID=ID, ref_dods = ref_dods, tidy=TRUE, tidy_style = "snake_case")
    suppressMessages(suppressWarnings(mnis_df <- dplyr::inner_join(mnis_df, opposition_posts_df)))}

    if (other_parliaments == TRUE) {
        other_parliaments_df <- mnis_other_parliaments(ID=ID, ref_dods = ref_dods, tidy=TRUE, tidy_style = "snake_case")
    suppressMessages(suppressWarnings(mnis_df <- dplyr::inner_join(mnis_df, other_parliaments_df)))}

    if (parliamentary_posts == TRUE) {
        parliamentary_posts_df <- mnis_parliamentary_posts(ID=ID, ref_dods = ref_dods, tidy=TRUE, tidy_style = "snake_case")
    suppressMessages(suppressWarnings(mnis_df <- dplyr::inner_join(mnis_df, parliamentary_posts_df)))}

    if (parties == TRUE) {
        parties_df <- mnis_parties(ID=ID, ref_dods = ref_dods, tidy=TRUE, tidy_style = "snake_case")
    suppressMessages(suppressWarnings(mnis_df <- dplyr::inner_join(mnis_df, parties_df)))}

    if (preferred_names == TRUE) {
        preferred_names_df <- mnis_preferred_names(ID=ID, ref_dods = ref_dods, tidy=TRUE, tidy_style = "snake_case")
    suppressMessages(suppressWarnings(mnis_df <- dplyr::inner_join(mnis_df, preferred_names_df)))}

    if (staff == TRUE) {
        staff_df <- mnis_staff(ID=ID, ref_dods = ref_dods, tidy=TRUE, tidy_style = "snake_case")
    suppressMessages(suppressWarnings(mnis_df <- dplyr::inner_join(mnis_df, staff_df)))}

    if (statuses == TRUE) {
        statuses_df <- mnis_statuses(ID=ID, ref_dods = ref_dods, tidy=TRUE, tidy_style = "snake_case")
    suppressMessages(suppressWarnings(mnis_df <- dplyr::inner_join(mnis_df, statuses_df)))}

    if (tidy == TRUE) {

      mnis_df <- mnis::mnis_tidy(mnis_df, tidy_style)

      mnis_df

    } else {

      mnis_df

    }

}

