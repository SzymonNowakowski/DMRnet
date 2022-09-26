release_questions <- function() {
  c(
    "Have you run the testing_branch tests?",
    "Have you updated DESCRIPTION with new version number?",
    "Have you updated README.md with new version number in master branch only (not testing_branch)?",
    "Have you updated NEWS.md with new version number?",
    "Have you updated .version.json with new version number?"
  )
}
